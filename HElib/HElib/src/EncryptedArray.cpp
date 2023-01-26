/* Copyright (C) 2012-2020 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */
/* EncryptedArray.cpp - Data-movement operations on arrays of slots
 */
#include <algorithm>
#include <helib/zzX.h>
#include <helib/EncryptedArray.h>
#include <helib/timing.h>
#include <helib/ClonedPtr.h>
#include <helib/norms.h>
#include <helib/exceptions.h>

#include "io.h"

namespace helib {

EncryptedArrayBase* buildEncryptedArray(const Context& context,
                                        const PAlgebraMod& alMod,
                                        const NTL::ZZX& G)
{
  if (alMod.getTag() == PA_cx_tag)
    return new EncryptedArrayCx(context, alMod.getCx());

  // By default
  // use the 1st factor F0
  const NTL::ZZX& GG = NTL::IsZero(G) ? alMod.getFactorsOverZZ()[0] : G;

  switch (alMod.getTag()) {
  case PA_GF2_tag: {
    return new EncryptedArrayDerived<PA_GF2>(context,
                                             NTL::conv<NTL::GF2X>(GG),
                                             alMod);
  }
  case PA_zz_p_tag: {
    NTL::zz_pBak bak;
    bak.save();
    alMod.restoreContext();
    return new EncryptedArrayDerived<PA_zz_p>(context,
                                              NTL::conv<NTL::zz_pX>(GG),
                                              alMod);
  }
  default:
    return nullptr;
  }
}

template <typename type>
EncryptedArrayDerived<type>::EncryptedArrayDerived(const Context& _context,
                                                   const RX& _G,
                                                   const PAlgebraMod& alMod) :
    context(_context), tab(alMod.getDerived(type()))
{
  tab.mapToSlots(mappingData, _G); // Compute the base-G representation maps
}

// rotate ciphertext in dimension i by amt
template <typename type>
void EncryptedArrayDerived<type>::rotate1D(Ctxt& ctxt,
                                           long i,
                                           long amt,
                                           bool dc) const
{
  HELIB_TIMER_START;
  helib::assertEq(&context, &ctxt.getContext(), "Context mismatch");
  helib::assertInRange(i,
                       0l,
                       dimension(),
                       "i must be between 0 and dimension()");

  RBak bak;
  bak.save();
  tab.restoreContext();

  const std::vector<std::vector<RX>>& maskTable = tab.getMaskTable();
  const PAlgebra& zMStar = getPAlgebra();
  long ord = sizeOfDimension(i);

  amt %= ord; // assumes division w/ remainder follows C++11
  if (amt == 0)
    return;
  if (amt < 0)
    amt += ord; // Make sure amt is in the range [1,ord-1]

  if (dc || nativeDimension(i)) { // native dimension or don't-care
    // For don't-care, we assume that any shifts "off the end" are zero
    ctxt.smartAutomorph(zMStar.genToPow(i, amt));
    return;
  }

  // more expensive "non-native" rotation

  helib::assertTrue(maskTable[i].size() > 0,
                    "Found non-positive sized mask table entry");

  ctxt.smartAutomorph(zMStar.genToPow(i, amt));
  // ctxt = \rho_i^{amt}(originalCtxt)

  Ctxt T(ctxt);
  T.smartAutomorph(zMStar.genToPow(i, -ord));
  // T = \rho_i^{amt-ord}(originalCtxt).
  // This strategy is geared toward the
  // assumption that we have the key switch matrix
  // for \rho_i^{-ord}

  const RX& mask = maskTable[i][amt];
  zzX mask_poly = balanced_zzX(mask);
  double sz = embeddingLargestCoeff(mask_poly, zMStar);
  DoubleCRT m1(mask_poly, context, ctxt.getPrimeSet() | T.getPrimeSet());
  // m1 will be used to multiply both ctxt and T

  // Compute ctxt = ctxt*m1 + T - T*m1
  ctxt.multByConstant(m1, sz);
  ctxt += T;
  T.multByConstant(m1, sz);
  ctxt -= T;
}

// Shift k positions along the i'th dimension with zero fill.
// Negative shift amount denotes shift in the opposite direction.
template <typename type>
void EncryptedArrayDerived<type>::shift1D(Ctxt& ctxt, long i, long k) const
{
  HELIB_TIMER_START;
  const PAlgebra& al = getPAlgebra();

  const std::vector<std::vector<RX>>& maskTable = tab.getMaskTable();

  RBak bak;
  bak.save();
  tab.restoreContext();

  assertEq(&context, &ctxt.getContext(), "Context mismatch");
  assertInRange(i,
                0l,
                (long)(al.numOfGens()),
                "i must be non-negative and less than the PAlgebra's "
                "generator count");

  long ord = al.OrderOf(i);

  if (k <= -ord || k >= ord) {
    ctxt.clear();
    return;
  }

  // Make sure amt is in the range [1,ord-1]
  long amt = k % ord;
  if (amt == 0)
    return;
  if (amt < 0)
    amt += ord;

  RX mask = maskTable[i][ord - amt];

  long val;
  if (k < 0)
    val = al.genToPow(i, amt - ord);
  else {
    mask = 1 - mask;
    val = al.genToPow(i, amt);
  }
  ctxt.multByConstant(balanced_zzX(mask)); // zero out slots where mask=0
  ctxt.smartAutomorph(val);                // shift left by val
  HELIB_TIMER_STOP;
}

// NOTE: masking depth: if there are N dimensions, and if for i = 1..N
// we define c_i = 1 if dimension i is bad and 0 o/w, then the masking
// depth is N - 1 + \sum_{i=1} c_i.

template <typename type>
void EncryptedArrayDerived<type>::rotate(Ctxt& ctxt, long amt) const
{
  HELIB_TIMER_START;

  const PAlgebra& al = getPAlgebra();
  const std::vector<std::vector<RX>>& maskTable = tab.getMaskTable();

  RBak bak;
  bak.save();
  tab.restoreContext();

  assertEq(&context, &ctxt.getContext(), "Context mismatch");

  // Simple case: just one generator
  if (al.numOfGens() == 1) { // VJS: bug fix: <= must be ==
    rotate1D(ctxt, 0, amt);
    return;
  }

  // Make sure that amt is in [1,nslots-1]
  amt %= (long)al.getNSlots();
  if (amt == 0) {
    return;
  }
  if (amt < 0)
    amt += al.getNSlots();

  // rotate the ciphertext, one dimension at a time
  long i = al.numOfGens() - 1;
  long v = al.coordinate(i, amt);
  RX mask = maskTable[i][v];
  Ctxt tmp(ctxt.getPubKey());
  const RXModulus& PhimXmod = tab.getPhimXMod();

  // optimize for the common case where the last generator has order in
  // Zm*/(p) different than its order in Zm*. In this case we can combine
  // the rotate1D relative to this generator with the masking after the
  // rotation. This saves one mult-by-constant, since we use the same mask
  // inside rotate1D as in the loop below.

  if (al.SameOrd(i) || v == 0)
    rotate1D(ctxt, i, v); // no need to optimize
  else {

    long ord = al.OrderOf(i);

    ctxt.smartAutomorph(al.genToPow(i, v));
    // ctxt = \rho_i^{v}(originalCtxt)

    tmp = ctxt;
    tmp.smartAutomorph(al.genToPow(i, -ord));
    // tmp = \rho_i^{v-ord}(originalCtxt).
    // This strategy assumes is geared toward the
    // assumption that we have the key switch matrix
    // for \rho_i^{-ord}

    zzX mask_poly = balanced_zzX(mask);
    double sz = embeddingLargestCoeff(mask_poly, al);

    DoubleCRT m1(mask_poly, context, ctxt.getPrimeSet() | tmp.getPrimeSet());
    // m1 will be used to multiply both ctxt and tmp

    // Compute ctxt = ctxt*m1, tmp = tmp*(1-m1)
    ctxt.multByConstant(m1, sz);

    Ctxt tmp1(tmp);
    tmp1.multByConstant(m1, sz);
    tmp -= tmp1;

    // apply rotation relative to next generator before combining the parts
    --i;
    v = al.coordinate(i, amt);
    rotate1D(ctxt, i, v);
    rotate1D(tmp, i, v + 1);
    ctxt += tmp; // combine the two parts

    if (i <= 0) {
      return;
    } // no more generators

    // update the mask for next iteration
    mask = ((mask * (maskTable[i][v] - maskTable[i][v + 1])) % PhimXmod) +
           maskTable[i][v + 1];
  }

  // Handle rotation relative to all the other generators (if any)
  for (i--; i >= 0; i--) {
    v = al.coordinate(i, amt);

    zzX mask_poly = balanced_zzX(mask);

    tmp = ctxt;
    tmp.multByConstant(mask_poly); // only the slots in which mask=1
    ctxt -= tmp;                   // only the slots in which mask=0

    rotate1D(tmp, i, v);
    rotate1D(ctxt, i, v + 1);
    ctxt += tmp;
    if (i > 0) {
      mask = ((mask * (maskTable[i][v] - maskTable[i][v + 1])) % PhimXmod) +
             maskTable[i][v + 1]; // update the mask for next iteration
    }
  }
  HELIB_TIMER_STOP;
}

template <typename type>
void EncryptedArrayDerived<type>::shift(Ctxt& ctxt, long k) const
{
  HELIB_TIMER_START;

  const PAlgebra& al = getPAlgebra();

  const std::vector<std::vector<RX>>& maskTable = tab.getMaskTable();

  RBak bak;
  bak.save();
  tab.restoreContext();

  assertEq(&context, &ctxt.getContext(), "Context mismatch");

  // Simple case: just one generator
  if (al.numOfGens() == 1) {
    shift1D(ctxt, 0, k);
    return;
  }

  long nSlots = al.getNSlots();

  // Shifting by more than the number of slots gives an all-zero ciphertext
  if (k <= -nSlots || k >= nSlots) {
    ctxt.multByConstant(NTL::to_ZZ(0));
    return;
  }

  // Make sure that amt is in [1,nslots-1]
  long amt = k % nSlots;
  if (amt == 0)
    return;
  if (amt < 0)
    amt += nSlots;

  // rotate the ciphertext, one dimension at a time
  long i = al.numOfGens() - 1;
  long v = al.coordinate(i, amt);
  RX mask = maskTable[i][v];
  Ctxt tmp(ctxt.getPubKey());
  const RXModulus& PhimXmod = tab.getPhimXMod();

  rotate1D(ctxt, i, v);
  for (i--; i >= 0; i--) {
    v = al.coordinate(i, amt);

    zzX mask_poly = balanced_zzX(mask);

    tmp = ctxt;
    tmp.multByConstant(mask_poly); // only the slots in which mask=1
    ctxt -= tmp;                   // only the slots in which mask=0
    if (i > 0) {
      rotate1D(ctxt, i, v + 1);
      rotate1D(tmp, i, v);
      ctxt += tmp; // combine the two parts

      mask = ((mask * (maskTable[i][v] - maskTable[i][v + 1])) % PhimXmod) +
             maskTable[i][v + 1]; // update the mask before next iteration
    } else {                      // i == 0
      if (k < 0)
        v -= al.OrderOf(0);
      shift1D(tmp, 0, v);
      shift1D(ctxt, 0, v + 1);
      ctxt += tmp;
    }
  }
  HELIB_TIMER_STOP;
}

template <typename type>
void EncryptedArrayDerived<type>::encode(NTL::ZZX& ptxt,
                                         const std::vector<RX>& array) const
{
  RX pp;
  tab.embedInSlots(pp, array, mappingData);

  // NOTE: previous version was
  //   ptxt = conv<NTL::ZZX>(pp);
  // which did not do balanced remainders at all
  zzX pp1 = balanced_zzX(pp);
  convert(ptxt, pp1);
}

template <typename type>
void EncryptedArrayDerived<type>::decode(std::vector<RX>& array,
                                         const NTL::ZZX& ptxt) const
{
  HELIB_TIMER_START;
  RX pp;
  conv(pp, ptxt);
  tab.decodePlaintext(array, pp, mappingData);
  HELIB_TIMER_STOP;
}

template <typename type>
void EncryptedArrayDerived<type>::encode(RX& ptxt,
                                         const std::vector<RX>& array) const
{
  tab.embedInSlots(ptxt, array, mappingData);
}

template <typename type>
void EncryptedArrayDerived<type>::decode(std::vector<RX>& array,
                                         const RX& ptxt) const
{
  tab.decodePlaintext(array, ptxt, mappingData);
}

template <typename type>
void EncryptedArrayDerived<type>::encode(NTL::ZZX& ptxt,
                                         const PlaintextArray& array) const
{
  RBak bak;
  bak.save();
  tab.restoreContext();
  encode(ptxt, array.getData<type>());
}

template <typename type>
void EncryptedArrayDerived<type>::decode(PlaintextArray& array,
                                         const NTL::ZZX& ptxt) const
{
  RBak bak;
  bak.save();
  tab.restoreContext();
  decode(array.getData<type>(), ptxt);
}

template <typename type>
void EncryptedArrayDerived<type>::encodeUnitSelector(zzX& ptxt, long i) const
{
  assertInRange(
      i,
      0l,
      (long)getPAlgebra().getNSlots(),
      "i must be non-negative and less than the PAlgebra's slot count");
  RBak bak;
  bak.save();
  tab.restoreContext();
  RX res;
  div(res, tab.getPhimXMod(), tab.getFactors()[i]);
  mul(res, res, tab.getCrtCoeffs()[i]);

  ptxt = balanced_zzX(res);
  // NOTE: previous version was
  //   convert(ptxt, res);
  // which did not do properly balanced remainders in some cases
}

template <typename type>
void EncryptedArrayDerived<type>::encode(zzX& ptxt,
                                         const std::vector<RX>& array) const
{
  RX pp;
  tab.embedInSlots(pp, array, mappingData);

  // NOTE: previous version was
  //   convert(ptxt, pp);
  // which did not do properly balanced remainders in some cases
  ptxt = balanced_zzX(pp);
}

template <typename type>
void EncryptedArrayDerived<type>::encode(zzX& ptxt,
                                         const PlaintextArray& array) const
{
  RBak bak;
  bak.save();
  tab.restoreContext();
  encode(ptxt, array.getData<type>());
}

template <typename type>
void EncryptedArrayDerived<type>::decode(std::vector<RX>& array,
                                         const NTL::Vec<long>& ptxt) const
{
  HELIB_TIMER_START;
  RX pp;
  convert(pp, ptxt);
  tab.decodePlaintext(array, pp, mappingData);
  HELIB_TIMER_STOP;
}

template <typename type>
void EncryptedArrayDerived<type>::decode(PlaintextArray& array,
                                         const NTL::Vec<long>& ptxt) const
{
  RBak bak;
  bak.save();
  tab.restoreContext();
  decode(array.getData<type>(), ptxt);
}

// this routine generates a "random" normal element and initializes a
// matrix mapping from polynomial to normal basis and its inverse. It
// uses a randomized algorithm to find the normal basis, but we init
// the PRG seed deterministically to ensure that we always get the
// same one (for a given set of parameters)

template <typename type>
void EncryptedArrayDerived<type>::initNormalBasisMatrix() const
{
  RandomState state;
  SetSeed(NTL::to_ZZ(1));
  do {
    typename NTL::Lazy<NTL::Pair<NTL::Mat<R>, NTL::Mat<R>>>::Builder builder(
        normalBasisMatrices);

    if (!builder())
      break;

    RBak bak;
    bak.save();
    restoreContext();
    REBak ebak;
    ebak.save();
    restoreContextForG();

    long d = RE::degree();
    long p = getPAlgebra().getP();
    long r = tab.getR();

    // compute change of basis matrix CB
    mat_R CB;
    CB.SetDims(d, d);
    RE normal_element;
    RE H;
    bool got_it = false;

    H = power(NTL::conv<RE>(RX(1, 1)), p);

    do {
      NTL::random(normal_element);

      RE pow;
      pow = normal_element;
      VectorCopy(CB[0], rep(pow), d);
      for (long i = 1; i < d; i++) {
        pow = eval(rep(pow), H);
        VectorCopy(CB[i], rep(pow), d);
      }

      NTL::Mat<NTL::ZZ> CB1;
      conv(CB1, CB);

      {
        NTL::zz_pBak bak1;
        bak1.save();
        NTL::zz_p::init(p);
        NTL::Mat<NTL::zz_p> CB2;
        conv(CB2, CB1);
        got_it = determinant(CB2) != 0;
      }
    } while (!got_it);

    NTL::Mat<R> CBi;
    ppInvert(CBi, CB, p, r);

    NTL::UniquePtr<NTL::Pair<NTL::Mat<R>, NTL::Mat<R>>> ptr;
    ptr.make(CB, CBi);
    builder.move(ptr);
  } while (0);
}

// PtxtArray member functions

void PtxtArray::writeToJSON(std::ostream& os) const
{
  executeRedirectJsonError<void>([&]() { os << this->writeToJSON(); });
}

JsonWrapper PtxtArray::writeToJSON() const
{
  auto body = [this]() {
    json jslots;

    if (ea.isCKKS()) {
      // When it is CKKS
      std::vector<std::complex<double>> data;
      store(data);
      jslots = data;
    } else {
      // When is it BGV
      std::vector<NTL::ZZX> data;
      store(data);
      std::vector<std::vector<long>> slots(data.size());
      for (std::size_t i = 0; i < data.size(); ++i) {
        long deg = NTL::deg(data[i]);
        if (deg == -1) {
          slots[i].emplace_back(0);
        }
        for (long j = 0; j <= deg; ++j) {
          slots[i].emplace_back(NTL::conv<long>(data[i][j]));
        }
      }
      jslots = slots;
    }

    json j{{"scheme", (ea.isCKKS() ? "CKKS" : "BGV")}, {"slots", jslots}};

    return wrap(toTypedJson<PtxtArray>(j));
  };

  return executeRedirectJsonError<JsonWrapper>(body);
}

PtxtArray PtxtArray::readFromJSON(std::istream& is, const Context& context)
{
  PtxtArray ret{context};
  ret.readJSON(is);
  return ret;
}

PtxtArray PtxtArray::readFromJSON(const JsonWrapper& tjw,
                                  const Context& context)
{
  PtxtArray ret{context};
  ret.readJSON(tjw);
  return ret;
}

void PtxtArray::readJSON(std::istream& is)
{
  executeRedirectJsonError<void>([&]() {
    json j;
    is >> j;
    this->readJSON(wrap(j));
  });
}

void PtxtArray::readJSON(const JsonWrapper& tjw)
{
  auto body = [&]() {
    json tj = unwrap(tjw);
    json jslots;
    // if the input is just an array short-circuit to slot deserialization
    // (assuming there is no type-header).
    if (tj.is_array()) {
      jslots = tj;

    } else {
      json j = fromTypedJson<PtxtArray>(tj);

      std::string expected_scheme{j.at("scheme").get<std::string>()};
      assertTrue<IOError>(
          expected_scheme == (ea.isCKKS() ? "CKKS" : "BGV"),
          "Scheme mismatch in deserialization.\nExpected: " + expected_scheme +
              ", actual: " + std::string(ea.isCKKS() ? "CKKS" : "BGV") + ".");

      jslots = j.at("slots");

      if (!jslots.is_array()) {
        throw IOError("Slot content is not a JSON array");
      }
    }

    if (static_cast<long>(jslots.size()) > this->getEA().size()) {
      std::stringstream err_msg;
      err_msg << "Cannot deserialize to PtxtArray: not enough slots.  "
              << "Trying to deserialize " << jslots.size() << " elements.  "
              << "Got " << this->getEA().size() << " slots.";
      throw IOError(err_msg.str());
    }

    if (ea.isCKKS()) {
      // Scheme is CKKS
      this->load(jslots.get<std::vector<std::complex<double>>>());
    } else {
      // Scheme is BGV
      auto json2data = [](const json& jslots) {
        std::vector<NTL::ZZX> data;
        data.reserve(jslots.size());
        for (const auto& jslot : jslots) {
          NTL::ZZX slot;
          if (jslot.is_array()) {
            for (std::size_t i = 0; i < jslot.size(); ++i) {
              NTL::SetCoeff(slot, i, static_cast<long>(jslot[i]));
            }
          } else {
            // Slot is a single number
            slot = static_cast<long>(jslot);
          }
          data.emplace_back(slot);
        }
        return data;
      };
      this->load(json2data(jslots));
    }
  };

  executeRedirectJsonError<void>(body);
}

std::istream& operator>>(std::istream& is, PtxtArray& pa)
{
  pa.readJSON(is);
  return is;
}

std::ostream& operator<<(std::ostream& os, const PtxtArray& pa)
{
  pa.writeToJSON(os);
  return os;
}

// Other functions...

void runningSums(const EncryptedArray& ea, Ctxt& ctxt)
{
  long n = ea.size();

  long shamt = 1;
  while (shamt < n) {
    Ctxt tmp = ctxt;
    ea.shift(tmp, shamt);
    ctxt += tmp; // ctxt = ctxt + (ctxt >> shamt)
    shamt = 2 * shamt;
  }
}

void totalSums(const EncryptedArray& ea, Ctxt& ctxt)
{
  long n = ea.size();

  if (n == 1)
    return;

  Ctxt orig = ctxt;

  long k = NTL::NumBits(n);
  long e = 1;

  for (long i = k - 2; i >= 0; i--) {
    Ctxt tmp1 = ctxt;
    ea.rotate(tmp1, e);
    ctxt += tmp1; // ctxt = ctxt + (ctxt >>> e)
    e = 2 * e;

    if (NTL::bit(n, i)) {
      Ctxt tmp2 = orig;
      ea.rotate(tmp2, e);
      ctxt += tmp2; // ctxt = ctxt + (orig >>> e)
                    // NOTE: we could have also computed
                    // ctxt =  (ctxt >>> e) + orig, however,
                    // this would give us greater depth/noise
      e += 1;
    }
  }
}

// Linearized polynomials.
// L describes a linear map M by describing its action on the standard
// power basis: M(x^j mod G) = (L[j] mod G), for j = 0..d-1.
// The result is a coefficient vector C for the linearized polynomial
// representing M: a polynomial h in Z/(p^r)[X] of degree < d is sent to
//
//    M(h(X) \bmod G)= \sum_{i=0}^{d-1}(C[j] \cdot h(X^{p^j}))\bmod G).
template <typename type>
void EncryptedArrayDerived<type>::buildLinPolyCoeffs(
    std::vector<NTL::ZZX>& C,
    const std::vector<NTL::ZZX>& L) const
{
  RBak bak;
  bak.save();
  restoreContext();
  std::vector<RX> CC, LL;
  convert(LL, L);
  buildLinPolyCoeffs(CC, LL);
  convert(C, CC);
}

template <typename type>
void EncryptedArrayDerived<type>::buildLinPolyCoeffs(
    std::vector<RX>& C,
    const std::vector<RX>& L) const
{
  HELIB_TIMER_START;

  RBak bak;
  bak.save();
  restoreContext(); // the NTL context for mod p^r
  REBak ebak;
  ebak.save();
  restoreContextForG(); // The NTL context for mod G

  do {
    typename NTL::Lazy<NTL::Mat<RE>>::Builder builder(linPolyMatrix);
    if (!builder())
      break;

    HELIB_NTIMER_START(buildLinPolyCoeffs_invert);

    long p = getPAlgebra().getP();
    long r = tab.getR();

    NTL::Mat<RE> M1;
    // build d x d matrix, d is taken from the current NTL context for G
    buildLinPolyMatrix(M1, p);
    NTL::Mat<RE> M2;
    ppInvert(M2, M1, p, r); // invert modulo prime-power p^r

    NTL::UniquePtr<NTL::Mat<RE>> ptr;
    ptr.make(M2);
    builder.move(ptr);
  } while (0);

  NTL::Vec<RE> CC, LL;
  convert(LL, L);
  mul(CC, LL, *linPolyMatrix);
  convert(C, CC);
}

// Apply the same linear transformation to all the slots.
// C[0...d-1] is the output of ea.buildLinPolyCoeffs
void applyLinPoly1(const EncryptedArray& ea,
                   Ctxt& ctxt,
                   const std::vector<NTL::ZZX>& C)
{
  assertEq(&ea.getContext(), &ctxt.getContext(), "Context mismatch");
  long d = ea.getDegree();
  assertEq(d, lsize(C), "ea's degree does not match the size of C");

  long nslots = ea.size();

  std::vector<NTL::ZZX> encodedC(d);
  for (long j = 0; j < d; j++) {
    std::vector<NTL::ZZX> v(nslots); // all the slots of v equal C[j]
    for (long i = 0; i < nslots; i++)
      v[i] = C[j];
    ea.encode(encodedC[j], v);
  }

  applyLinPolyLL(ctxt, encodedC, ea.getDegree());
}

// Apply different transformations to different slots. Each row in
// the matrix Cvec[0...nslots-1][0...d-1] is a length-d vector which
// is the output of ea.buildLinPolyCoeffs
void applyLinPolyMany(const EncryptedArray& ea,
                      Ctxt& ctxt,
                      const std::vector<std::vector<NTL::ZZX>>& Cvec)
{
  assertEq(&ea.getContext(), &ctxt.getContext(), "Context mismatch");
  long d = ea.getDegree();
  long nslots = ea.size();

  assertEq(nslots, lsize(Cvec), "Number of slots does not match size of Cvec");
  for (long i = 0; i < nslots; i++) {
    assertEq(d,
             lsize(Cvec[i]),
             "Found entry of Cvec with size unequal to degree of ea");
  }

  std::vector<NTL::ZZX> encodedC(d);
  for (long j = 0; j < d; j++) {     // encodedC[j] encodes j'th column in Cvec
    std::vector<NTL::ZZX> v(nslots); // copy j'th column to v
    for (long i = 0; i < nslots; i++)
      v[i] = Cvec[i][j];
    ea.encode(encodedC[j], v); // then encode it
  }

  applyLinPolyLL(ctxt, encodedC, ea.getDegree());
}

// A low-level variant: encodedCoeffs has all the linPoly coeffs encoded
// in slots; different transformations can be encoded in different slots
template <typename P>
void applyLinPolyLL(Ctxt& ctxt, const std::vector<P>& encodedC, long d)
{
  assertEq(d, lsize(encodedC), "d does not match size of encodedC");

  ctxt.cleanUp(); // not sure, but this may be a good idea

  Ctxt tmp(ctxt);

  ctxt.multByConstant(encodedC[0]);
  for (long j = 1; j < d; j++) {
    Ctxt tmp1(tmp);
    tmp1.frobeniusAutomorph(j);
    tmp1.multByConstant(encodedC[j]);
    ctxt += tmp1;
  }
}
template void applyLinPolyLL(Ctxt& ctxt,
                             const std::vector<zzX>& encodedC,
                             long d);
template void applyLinPolyLL(Ctxt& ctxt,
                             const std::vector<NTL::ZZX>& encodedC,
                             long d);
template void applyLinPolyLL(Ctxt& ctxt,
                             const std::vector<DoubleCRT>& encodedC,
                             long d);

/****************** End linear transformation code ******************/
/********************************************************************/

// PlaintextArray

template <typename type>
class rotate_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    long k)
  {
    PA_BOILER(type)

    std::vector<RX> tmp(n);

    for (long i = 0; i < n; i++)
      tmp[((i + k) % n + n) % n] = data[i];

    data = tmp;
  }
};

void rotate(const EncryptedArray& ea, PlaintextArray& pa, long k)
{
  ea.dispatch<rotate_pa_impl>(pa, k);
}

//=============================================================================

template <typename type>
class rotate1D_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    long i,
                    long k)
  {
    PA_BOILER(type)

    helib::assertInRange(i,
                         0l,
                         ea.dimension(),
                         "i must be between 0 and dimension()");

    std::vector<RX> tmp(n);
    ea.EncryptedArrayBase::rotate1D(tmp, data, i, k);
    data = tmp;
  }
};

void rotate1D(const EncryptedArray& ea, PlaintextArray& pa, long i, long k)
{
  ea.dispatch<rotate1D_pa_impl>(pa, i, k);
}

//=============================================================================

template <typename type>
class shift_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    long k)
  {
    PA_BOILER(type)

    for (long j : range(n))
      if (j + k >= n || j + k < 0)
        data[j] = 0;

    rotate_pa_impl<type>::apply(ea, pa, k);
  }
};

void shift(const EncryptedArray& ea, PlaintextArray& pa, long k)
{
  ea.dispatch<shift_pa_impl>(pa, k);
}

//=============================================================================

template <typename type>
class shift1D_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    long i,
                    long k)
  {
    PA_BOILER(type)

    helib::assertInRange(i,
                         0l,
                         ea.dimension(),
                         "i must be between 0 and dimension()");

    long sz = ea.sizeOfDimension(i);

    for (long j : range(n)) {
      long c = ea.coordinate(i, j);
      if (c + k >= sz || c + k < 0)
        data[j] = 0;
    }

    rotate1D_pa_impl<type>::apply(ea, pa, i, k);
  }
};

void shift1D(const EncryptedArray& ea, PlaintextArray& pa, long i, long k)
{
  ea.dispatch<shift1D_pa_impl>(pa, i, k);
}

//=============================================================================

template <typename type>
class encode_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    const std::vector<long>& array)
  {
    PA_BOILER(type)

    long m = std::min(n, lsize(array));
    for (long i = 0; i < m; i++)
      data[i] = array[i];
    for (long i = m; i < n; i++)
      data[i] = 0;
  }

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    const std::vector<NTL::ZZX>& array)
  {
    PA_BOILER(type)

    long m = std::min(n, lsize(array));
    for (long i = 0; i < m; i++) {
      if (deg(array[i]) < d) {
        convert(data[i], array[i]);
      } else {
        RX a;
        convert(a, array[i]);
        data[i] = a % G;
      }
    }
    for (long i = m; i < n; i++)
      data[i] = 0;
  }

  static void apply(UNUSED const EncryptedArrayDerived<type>& ea,
                    UNUSED PlaintextArray& pa,
                    UNUSED const std::vector<cx_double>& array)
  {
    throw LogicError("function not implemented");
  }

  static void apply(UNUSED const EncryptedArrayDerived<type>& ea,
                    UNUSED PlaintextArray& pa,
                    UNUSED const std::vector<double>& array)
  {
    throw LogicError("function not implemented");
  }
};

template <>
class encode_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    PlaintextArray& pa,
                    const std::vector<long>& array)
  {
    PA_BOILER(PA_cx)

    long m = std::min(n, lsize(array));
    for (long i = 0; i < m; i++)
      data[i] = array[i];
    for (long i = m; i < n; i++)
      data[i] = 0;
  }

  static void apply(UNUSED const EncryptedArrayDerived<PA_cx>& ea,
                    UNUSED PlaintextArray& pa,
                    UNUSED const std::vector<NTL::ZZX>& array)
  {
    throw LogicError("function not implemented");
  }

  static void apply(UNUSED const EncryptedArrayDerived<PA_cx>& ea,
                    UNUSED PlaintextArray& pa,
                    UNUSED const std::vector<cx_double>& array)
  {
    PA_BOILER(PA_cx)

    long m = std::min(n, lsize(array));
    for (long i = 0; i < m; i++)
      data[i] = array[i];
    for (long i = m; i < n; i++)
      data[i] = 0;
  }

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    PlaintextArray& pa,
                    const std::vector<double>& array)
  {
    PA_BOILER(PA_cx)

    long m = std::min(n, lsize(array));
    for (long i = 0; i < m; i++)
      data[i] = array[i];
    for (long i = m; i < n; i++)
      data[i] = 0;
  }
};

void encode(const EncryptedArray& ea,
            PlaintextArray& pa,
            const std::vector<long>& array)
{
  ea.dispatch<encode_pa_impl>(pa, array);
}

void encode(const EncryptedArray& ea,
            PlaintextArray& pa,
            const std::vector<NTL::ZZX>& array)
{
  ea.dispatch<encode_pa_impl>(pa, array);
}

void encode(const EncryptedArray& ea,
            PlaintextArray& pa,
            const std::vector<cx_double>& array)
{
  ea.dispatch<encode_pa_impl>(pa, array);
}

void encode(const EncryptedArray& ea,
            PlaintextArray& pa,
            const std::vector<double>& array)
{
  ea.dispatch<encode_pa_impl>(pa, array);
}

void encode(const EncryptedArray& ea, PlaintextArray& pa, long val)
{
  long n = ea.size();
  std::vector<long> array;
  array.resize(n);
  for (long i = 0; i < n; i++)
    array[i] = val;
  encode(ea, pa, array);
}

void encode(const EncryptedArray& ea, PlaintextArray& pa, const NTL::ZZX& val)
{
  long n = ea.size();
  std::vector<NTL::ZZX> array;
  array.resize(n);
  for (long i = 0; i < n; i++)
    array[i] = val;
  encode(ea, pa, array);
}

void encode(const EncryptedArray& ea, PlaintextArray& pa, cx_double val)
{
  long n = ea.size();
  std::vector<cx_double> array;
  array.resize(n);
  for (long i = 0; i < n; i++)
    array[i] = val;
  encode(ea, pa, array);
}

void encode(const EncryptedArray& ea, PlaintextArray& pa, double val)
{
  long n = ea.size();
  std::vector<cx_double> array;
  array.resize(n);
  for (long i = 0; i < n; i++)
    array[i] = val;
  encode(ea, pa, array);
}

//=============================================================================

template <typename type>
class randomReal_pa_impl
{
public:
  // TODO: The following may be removed as its injected fields are never used
  // PA_INJECT(type)

  static void apply(UNUSED const EncryptedArrayDerived<type>& ea,
                    UNUSED PlaintextArray& pa)
  {
    // TODO: The following may be removed as the function is not implemented
    // PA_BOILER(type)

    throw LogicError("function not implemented");
  }
};

template <>
class randomReal_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  static void apply(const EncryptedArrayDerived<PA_cx>& ea, PlaintextArray& pa)
  {
    PA_BOILER(PA_cx)

    for (long i = 0; i < n; i++)
      data[i] = RandomReal();
  }
};

void randomReal(const EncryptedArray& ea, PlaintextArray& pa)
{
  ea.dispatch<randomReal_pa_impl>(pa);
}

//=============================================================================

template <typename type>
class random_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, PlaintextArray& pa)
  {
    PA_BOILER(type)

    for (long i = 0; i < n; i++)
      random(data[i], d);
  }
};

template <>
class random_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  // same as randomReal
  static void apply(const EncryptedArrayDerived<PA_cx>& ea, PlaintextArray& pa)
  {
    PA_BOILER(PA_cx)

    for (long i = 0; i < n; i++)
      data[i] = RandomReal();
  }
};

void random(const EncryptedArray& ea, PlaintextArray& pa)
{
  ea.dispatch<random_pa_impl>(pa);
}

//=============================================================================

template <typename type>
class randomComplex_pa_impl
{
public:
  // TODO: The following may be removed as its injected fields are never used
  // PA_INJECT(type)

  static void apply(UNUSED const EncryptedArrayDerived<type>& ea,
                    UNUSED PlaintextArray& pa)
  {
    // TODO: The following may be removed as the function is not implemented
    // PA_BOILER(type)

    throw LogicError("function not implemented");
  }
};

template <>
class randomComplex_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  static void apply(const EncryptedArrayDerived<PA_cx>& ea, PlaintextArray& pa)
  {
    PA_BOILER(PA_cx)

    for (long i = 0; i < n; i++)
      data[i] = RandomComplex();
  }
};

void randomComplex(const EncryptedArray& ea, PlaintextArray& pa)
{
  ea.dispatch<randomComplex_pa_impl>(pa);
}

//=============================================================================

template <typename type>
class decode_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    std::vector<NTL::ZZX>& array,
                    const PlaintextArray& pa)
  {
    CPA_BOILER(type)

    convert(array, data);
  }

  static void apply(const EncryptedArrayDerived<type>& ea,
                    std::vector<long>& array,
                    const PlaintextArray& pa)
  {
    CPA_BOILER(type)

    convert(array, data);
  }

  static void apply(UNUSED const EncryptedArrayDerived<type>& ea,
                    UNUSED std::vector<cx_double>& array,
                    UNUSED const PlaintextArray& pa)
  {
    throw LogicError("function not implemented");
  }

  static void apply(UNUSED const EncryptedArrayDerived<type>& ea,
                    UNUSED std::vector<double>& array,
                    UNUSED const PlaintextArray& pa)
  {
    throw LogicError("function not implemented");
  }
};

template <>
class decode_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  static void apply(UNUSED const EncryptedArrayDerived<PA_cx>& ea,
                    UNUSED std::vector<NTL::ZZX>& array,
                    UNUSED const PlaintextArray& pa)
  {
    throw LogicError("function not implemented");
  }

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    std::vector<long>& array,
                    const PlaintextArray& pa)
  {
    CPA_BOILER(PA_cx)

    project_and_round(array, data);
  }

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    std::vector<cx_double>& array,
                    const PlaintextArray& pa)
  {
    CPA_BOILER(PA_cx)

    convert(array, data);
  }

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    std::vector<double>& array,
                    const PlaintextArray& pa)
  {
    CPA_BOILER(PA_cx)

    project(array, data);
  }
};

void decode(const EncryptedArray& ea,
            std::vector<long>& array,
            const PlaintextArray& pa)
{
  ea.dispatch<decode_pa_impl>(array, pa);
}

void decode(const EncryptedArray& ea,
            std::vector<NTL::ZZX>& array,
            const PlaintextArray& pa)
{
  ea.dispatch<decode_pa_impl>(array, pa);
}

void decode(const EncryptedArray& ea,
            std::vector<cx_double>& array,
            const PlaintextArray& pa)
{
  ea.dispatch<decode_pa_impl>(array, pa);
}

void decode(const EncryptedArray& ea,
            std::vector<double>& array,
            const PlaintextArray& pa)
{
  ea.dispatch<decode_pa_impl>(array, pa);
}

//=============================================================================

template <typename type>
class equals_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    bool& res,
                    const PlaintextArray& pa,
                    const PlaintextArray& other)
  {
    CPA_BOILER(type)

    const std::vector<RX>& odata = other.getData<type>();
    res = (data == odata);
  }

  static void apply(const EncryptedArrayDerived<type>& ea,
                    bool& res,
                    const PlaintextArray& pa,
                    const std::vector<long>& other)
  {
    CPA_BOILER(type)

    std::vector<RX> odata;
    convert(odata, other);
    res = (data == odata);
  }

  static void apply(const EncryptedArrayDerived<type>& ea,
                    bool& res,
                    const PlaintextArray& pa,
                    const std::vector<NTL::ZZX>& other)
  {
    CPA_BOILER(type)

    std::vector<RX> odata;
    convert(odata, other);
    res = (data == odata);
  }
};

template <>
class equals_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  // exact comparison...not very useful
  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    bool& res,
                    const PlaintextArray& pa,
                    const PlaintextArray& other)
  {
    CPA_BOILER(PA_cx)

    const std::vector<RX>& odata = other.getData<PA_cx>();
    res = (data == odata);
  }

  static void apply(UNUSED const EncryptedArrayDerived<PA_cx>& ea,
                    UNUSED bool& res,
                    UNUSED const PlaintextArray& pa,
                    UNUSED const std::vector<long>& other)
  {
    throw LogicError("function not implemented");
  }

  static void apply(UNUSED const EncryptedArrayDerived<PA_cx>& ea,
                    UNUSED bool& res,
                    UNUSED const PlaintextArray& pa,
                    UNUSED const std::vector<NTL::ZZX>& other)
  {
    throw LogicError("function not implemented");
  }
};

bool equals(const EncryptedArray& ea,
            const PlaintextArray& pa,
            const PlaintextArray& other)
{
  bool res;
  ea.dispatch<equals_pa_impl>(res, pa, other);
  return res;
}

bool equals(const EncryptedArray& ea,
            const PlaintextArray& pa,
            const std::vector<long>& other)
{
  bool res;
  ea.dispatch<equals_pa_impl>(res, pa, other);
  return res;
}

bool equals(const EncryptedArray& ea,
            const PlaintextArray& pa,
            const std::vector<NTL::ZZX>& other)
{
  bool res;
  ea.dispatch<equals_pa_impl>(res, pa, other);
  return res;
}

//=============================================================================

template <typename type>
class add_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    const PlaintextArray& other)
  {
    PA_BOILER(type)

    const std::vector<RX>& odata = other.getData<type>();

    for (long i = 0; i < n; i++)
      data[i] += odata[i];
  }
};

void add(const EncryptedArray& ea,
         PlaintextArray& pa,
         const PlaintextArray& other)
{
  ea.dispatch<add_pa_impl>(pa, other);
}

//=============================================================================

template <typename type>
class sub_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    const PlaintextArray& other)
  {
    PA_BOILER(type)

    const std::vector<RX>& odata = other.getData<type>();

    for (long i = 0; i < n; i++)
      data[i] -= odata[i];
  }
};

void sub(const EncryptedArray& ea,
         PlaintextArray& pa,
         const PlaintextArray& other)
{
  ea.dispatch<sub_pa_impl>(pa, other);
}

//=============================================================================

template <typename type>
class mul_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    const PlaintextArray& other)
  {
    PA_BOILER(type)

    const std::vector<RX>& odata = other.getData<type>();

    for (long i = 0; i < n; i++)
      data[i] = (data[i] * odata[i]) % G;
  }
};

template <>
class mul_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    PlaintextArray& pa,
                    const PlaintextArray& other)
  {
    PA_BOILER(PA_cx)

    const std::vector<RX>& odata = other.getData<PA_cx>();

    for (long i = 0; i < n; i++)
      data[i] = data[i] * odata[i];
  }
};

void mul(const EncryptedArray& ea,
         PlaintextArray& pa,
         const PlaintextArray& other)
{
  ea.dispatch<mul_pa_impl>(pa, other);
}

//=============================================================================

template <typename type>
class negate_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, PlaintextArray& pa)
  {
    PA_BOILER(type)

    for (long i = 0; i < n; i++)
      data[i] = -data[i];
  }
};

void negate(const EncryptedArray& ea, PlaintextArray& pa)
{
  ea.dispatch<negate_pa_impl>(pa);
}

//=============================================================================

template <typename type>
class frobeniusAutomorph_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    long j)
  {
    PA_BOILER(type)

    long p = ea.getPAlgebra().getP();

    j = mcMod(j, d);
    RX H = NTL::PowerMod(RX(1, 1), NTL::power_ZZ(p, j), G);

    for (long i = 0; i < n; i++)
      data[i] = NTL::CompMod(data[i], H, G);
  }

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    const NTL::Vec<long>& vec)
  {
    PA_BOILER(type)

    assertEq(vec.length(), n, "vec has incorrect length");

    long p = ea.getPAlgebra().getP();

    for (long i = 0; i < n; i++) {
      long j = mcMod(vec[i], d);
      RX H = NTL::PowerMod(RX(1, 1), NTL::power_ZZ(p, j), G);
      data[i] = NTL::CompMod(data[i], H, G);
    }
  }
};

template <>
class frobeniusAutomorph_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    PlaintextArray& pa,
                    long j)
  {
    PA_BOILER(PA_cx)

    if (j % 2 != 0)
      for (long i = 0; i < n; i++)
        data[i] = std::conj(data[i]);
  }

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    PlaintextArray& pa,
                    const NTL::Vec<long>& vec)
  {
    PA_BOILER(PA_cx)

    assertEq(vec.length(), n, "vec has incorrect length");

    for (long i = 0; i < n; i++) {
      if (vec[i] % 2 != 0)
        data[i] = std::conj(data[i]);
    }
  }
};

void frobeniusAutomorph(const EncryptedArray& ea, PlaintextArray& pa, long j)
{
  ea.dispatch<frobeniusAutomorph_pa_impl>(pa, j);
}

void frobeniusAutomorph(const EncryptedArray& ea,
                        PlaintextArray& pa,
                        const NTL::Vec<long>& vec)
{
  ea.dispatch<frobeniusAutomorph_pa_impl>(pa, vec);
}

//=============================================================================

template <typename type>
class extractRealPart_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const UNUSED EncryptedArrayDerived<type>& ea,
                    UNUSED PlaintextArray& pa)
  {
    throw LogicError("function not implemented");
  }
};

template <>
class extractRealPart_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  static void apply(const EncryptedArrayDerived<PA_cx>& ea, PlaintextArray& pa)
  {
    PA_BOILER(PA_cx)

    for (long i = 0; i < n; i++)
      data[i] = data[i].real();
  }
};

void extractRealPart(const EncryptedArray& ea, PlaintextArray& pa)
{
  ea.dispatch<extractRealPart_pa_impl>(pa);
}

//=============================================================================

template <typename type>
class extractImPart_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const UNUSED EncryptedArrayDerived<type>& ea,
                    UNUSED PlaintextArray& pa)
  {
    throw LogicError("function not implemented");
  }
};

template <>
class extractImPart_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  static void apply(const EncryptedArrayDerived<PA_cx>& ea, PlaintextArray& pa)
  {
    PA_BOILER(PA_cx)

    for (long i = 0; i < n; i++)
      data[i] = data[i].imag();
  }
};

void extractImPart(const EncryptedArray& ea, PlaintextArray& pa)
{
  ea.dispatch<extractImPart_pa_impl>(pa);
}

//=============================================================================

void power(const EncryptedArray& ea, PlaintextArray& pa, long e)
{
  if (e <= 1)
    return;

  PlaintextArray pwr = pa; // holds x^{2^i} in i+1'st iteration
  encode(ea, pa, 1L);      // set pa =1 in every slot
  while (e > 0) {
    if (e & 1)
      mul(ea, pa, pwr); // multiply if needed
    mul(ea, pwr, pwr);  // square
    e >>= 1;
  }
}

//=============================================================================

template <typename type>
class applyPerm_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    PlaintextArray& pa,
                    const NTL::Vec<long>& pi)
  {
    PA_BOILER(type)

    assertEq(pi.length(), n, "pi has incorrect length");

    std::vector<RX> tmp;
    tmp.resize(n);
    for (long i = 0; i < n; i++)
      tmp[i] = data[pi[i]];

    data = tmp;
  }
};

void applyPerm(const EncryptedArray& ea,
               PlaintextArray& pa,
               const NTL::Vec<long>& pi)
{
  ea.dispatch<applyPerm_pa_impl>(pa, pi);
}

//=============================================================================

template <typename type>
class print_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    std::ostream& s,
                    const PlaintextArray& pa)
  {
    CPA_BOILER(type)

    s << data;
  }
};

void print(const EncryptedArray& ea, std::ostream& s, const PlaintextArray& pa)
{
  ea.dispatch<print_pa_impl>(s, pa);
}

//=============================================================================

template <typename type>
class Norm_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea,
                    double& res,
                    const PlaintextArray& pa)
  {
    CPA_BOILER(type)

    res = 0;
    for (long i = 0; i < n; i++)
      if (data[i] != 0) {
        res = 1;
        return;
      }
  }

  static void apply(const EncryptedArrayDerived<type>& ea,
                    double& res,
                    const PlaintextArray& pa,
                    const PlaintextArray& other)
  {
    CPA_BOILER(type)

    const std::vector<RX>& odata = other.getData<type>();

    res = 0;
    for (long i = 0; i < n; i++)
      if (data[i] != odata[i]) {
        res = 1;
        return;
      }
  }
};

template <>
class Norm_pa_impl<PA_cx>
{
public:
  PA_INJECT(PA_cx)

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    double& res,
                    const PlaintextArray& pa)
  {
    CPA_BOILER(PA_cx)

    res = Norm(data);
  }

  static void apply(const EncryptedArrayDerived<PA_cx>& ea,
                    double& res,
                    const PlaintextArray& pa,
                    const PlaintextArray& other)
  {
    CPA_BOILER(PA_cx)

    const std::vector<RX>& odata = other.getData<PA_cx>();

    res = Distance(data, odata);
  }
};

double Norm(const EncryptedArray& ea, const PlaintextArray& pa)
{
  double res;
  ea.dispatch<Norm_pa_impl>(res, pa);
  return res;
}

double Distance(const EncryptedArray& ea,
                const PlaintextArray& pa,
                const PlaintextArray& other)
{
  double res;
  ea.dispatch<Norm_pa_impl>(res, pa, other);
  return res;
}

//=============================================================================

template <typename type>
class totalSums_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, PlaintextArray& pa)
  {
    PA_BOILER(type)

    RX sum(0);

    for (long i = 0; i < n; i++)
      sum += data[i];

    for (long i = 0; i < n; i++)
      data[i] = sum;
  }
};

void totalSums(const EncryptedArray& ea, PlaintextArray& pa)
{
  ea.dispatch<totalSums_pa_impl>(pa);
}

//=============================================================================

template <typename type>
class runningSums_pa_impl
{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, PlaintextArray& pa)
  {
    PA_BOILER(type)

    for (long i = 1; i < n; i++)
      data[i] += data[i - 1];
  }
};

void runningSums(const EncryptedArray& ea, PlaintextArray& pa)
{
  ea.dispatch<runningSums_pa_impl>(pa);
}

//=============================================================================

// Explicit instantiation

template class EncryptedArrayDerived<PA_GF2>;
template class EncryptedArrayDerived<PA_zz_p>;

template class PlaintextArrayDerived<PA_GF2>;
template class PlaintextArrayDerived<PA_zz_p>;
template class PlaintextArrayDerived<PA_cx>;

} // namespace helib
