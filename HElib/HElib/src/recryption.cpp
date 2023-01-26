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
#include <NTL/BasicThreadPool.h>

#include <helib/recryption.h>
#include <helib/EncryptedArray.h>
#include <helib/EvalMap.h>
#include <helib/powerful.h>
#include <helib/CtPtrs.h>
#include <helib/intraSlot.h>
#include <helib/norms.h>
#include <helib/sample.h>
#include <helib/debugging.h>
#include <helib/fhe_stats.h>
#include <helib/log.h>
#include <helib/polyEval.h>

#include <algorithm>
#include <math.h>
#include <cmath>
#include <chrono>

#ifdef HELIB_DEBUG

#include <helib/debugging.h>

namespace helib {

long printFlag = FLAG_PRINT_VEC;
/************************ Some local functions ***********************/
/*********************************************************************/

static void checkCriticalValue(const std::vector<NTL::ZZX>& zzParts,
                               const DoubleCRT& sKey,
                               const RecryptData& rcData,
                               long q);

static void checkRecryptBounds(const std::vector<NTL::ZZX>& zzParts,
                               const DoubleCRT& sKey,
                               const Context& context,
                               long q);

static void checkRecryptBounds_v(const std::vector<NTL::ZZX>& v,
                                 const DoubleCRT& sKey,
                                 const Context& context,
                                 long q);
} // namespace helib

#endif // HELIB_DEBUG

namespace helib {

#include "Magma_declare.h"

// Return in poly a polynomial with X^i encoded in all the slots
static void x2iInSlots(NTL::ZZX& poly,
                       long i,
                       std::vector<NTL::ZZX>& xVec,
                       const EncryptedArray& ea)
{
  xVec.resize(ea.size());
  NTL::ZZX x2i = NTL::ZZX(i, 1);
  for (long j = 0; j < (long)xVec.size(); j++)
    xVec[j] = x2i;
  ea.encode(poly, xVec);
}

// Make every entry of vec divisible by p2e by adding/subtracting q, while
// keeping the added multiples small.  Specifically, for q = 1 mod p2e and any
// integer z can be made divisible by p2e via z' = z + v*q, with |v| <= p2e/2.

static void newMakeDivisible(NTL::ZZX& poly,
                             long p2e,
                             long q,
                             const Context& context,
                             NTL::ZZX& vpoly)
{
  if (p2e == 1) {
    vpoly = 0;
    return;
  }

  assertTrue<InvalidArgument>(q > 0l, "q must be positive");
  assertTrue<InvalidArgument>(p2e > 0l, "p2e must be positive");

  assertEq<InvalidArgument>(q % p2e, 1l, "q must equal 1 modulo p2e");

  long p = context.getP();

  const RecryptData& rcData = context.getRcData();
  const PowerfulDCRT& p2d_conv = *rcData.p2dConv;

  NTL::Vec<NTL::ZZ> pwrfl;
  p2d_conv.ZZXtoPowerful(pwrfl, poly);

#ifdef HELIB_DEBUG
  NTL::Vec<NTL::ZZ> vvec(NTL::INIT_SIZE, pwrfl.length());
#endif

  for (long i : range(pwrfl.length())) {
    NTL::ZZ& z = pwrfl[i];
    long v;

    // What to add to z to make it divisible by p2e?
    long zMod = rem(z, p2e); // zMod is in [0,p2e-1]
    // NOTE: this makes sure we get a truly balanced remainder
    if (zMod > p2e / 2 || (p == 2 && zMod == p2e / 2 && NTL::RandomBnd(2))) {
      // randomize so that v has expected value 0
      zMod = p2e - zMod;
    } else {
      // need to add a negative number
      zMod = -zMod;
    }
    v = zMod;
    z += NTL::to_ZZ(q) * v; // make z divisible by p2e

    if (rem(z, p2e) != 0) { // sanity check
      std::cerr << "**error: original z[" << i
                << "]=" << (z - (NTL::to_ZZ(q) * v)) << std::dec
                << ", p^e=" << p2e << std::endl;
      std::cerr << "z' = z + " << v << "*q = " << z << std::endl;
      exit(1);
    }

#ifdef HELIB_DEBUG
    vvec[i] = v;
#endif
  }

  p2d_conv.powerfulToZZX(poly, pwrfl);

#ifdef HELIB_DEBUG
  p2d_conv.powerfulToZZX(vpoly, vvec);
#endif
}

/*********************************************************************/
/*********************************************************************/

/**
 * Summary of Appendix A from https://ia.cr/2014/873 (version from 2019):
 * Assume that we already chosen e, e' and t.
 *
 * Based in this analysis, we need
 *    (1) (f*p^{e'} + 2*p^r+2))*B <= p^e/2
 * where B is a certain high-probability bound and f is a certain
 * fudge factor.
 *
 **/

// the routine compute_fudge is used to correct for the fact that
// the v-coeffs are not quite uniform

static double compute_fudge(long p2ePrime, long p2e)
{
  double eps = 0;

  if (p2ePrime > 1) {

    if (p2ePrime % 2 == 0) {
      eps = 1 / fsquare(p2ePrime);

      // The exact variance in this case is at most the variance
      // of a random variable that is distributed over
      //    -N..+N
      // where N = 2^{e'}/2.
      // Each endpoint occurs with probability 1/(4*N),
      // and the remaining values each occur with the same probability
      // 1/(2*N)

      // This variance is exactly computed as
      //    (N^2)/3 + 1/6 = ((N^2)/3)*(1 + 1/(2*N^2)), where N = 2^{e'}/2
      // So the std dev is at most
      //    N/sqrt(3)*(1 + 1/(4*N^2))

    } else {
      eps = 1 / double(p2e);

      // We are computing X + Y mod p^{e'}, where
      // X and Y are independent.
      // Y is uniformly distributed over
      //    -floor(p^{r}/2)..floor(p^{r}/2)
      // X is distributed over
      //    -floor(p^e/2)-1..floor(p^e/2)+1,
      // where each endpoint occurs with probability 1 / (2*(p^e+1)),
      // and the remaining p^e values are equally likely

      // The variance in this case is bounded by
      //   (N^2)/3*(1-eps) + (N^2)*eps = (N^2)/3*(1+2*eps),
      //       where N = p^{e'}/2 and eps < 1/p^e
      // So the std dev is bounded by
      //    N/sqrt(3)*sqrt(1+2*eps) <= N/sqrt(3)*(1+eps)
    }
  }

  return 1 + eps;
}

void RecryptData::setAE(long& e, long& ePrime, const Context& context)
{
  double coeff_bound = context.boundForRecryption();
  // coeff_bound is ultimately a high prob bound on |w0+w1*s|,
  // the coeffs of w0, w1 are chosen uniformly on [-1/2,1/2]

  long p = context.getP();
  long p2r = context.getAlMod().getPPowR();
  long r = context.getAlMod().getR();
  long frstTerm = 2 * p2r + 2;

  long e_bnd = 0;
  long p2e_bnd = 1;
  while (p2e_bnd <= ((1L << 30) - 2) / p) { // NOTE: this avoids overflow
    e_bnd++;
    p2e_bnd *= p;
  }
  // e_bnd is largest e such that p^e+1 < 2^30

  // Start with the smallest e s.t. p^e/2 >= frstTerm*coeff_bound
  ePrime = 0;
  e = r + 1;
  while (e <= e_bnd && NTL::power_long(p, e) < frstTerm * coeff_bound * 2)
    e++;

  //  if (e > e_bnd) Error("setAE: cannot find suitable e");
  assertFalse<RuntimeError>(e > e_bnd, "setAE: cannot find suitable e");

  // long ePrimeTry = r+1;
  long ePrimeTry = 1;

  while (ePrimeTry <= e_bnd) {
    long p2ePrimeTry = NTL::power_long(p, ePrimeTry);
    // long eTry = ePrimeTry+1;
    long eTry = std::max(r + 1, ePrimeTry + 1);
    while (eTry <= e_bnd && eTry - ePrimeTry < e - ePrime) {
      long p2eTry = NTL::power_long(p, eTry);
      double fudge = compute_fudge(p2ePrimeTry, p2eTry);
      if (p2eTry >= (p2ePrimeTry * fudge + frstTerm) * coeff_bound * 2)
        break;

      eTry++;
    }

    if (eTry <= e_bnd && eTry - ePrimeTry < e - ePrime) {
      e = eTry;
      ePrime = ePrimeTry;
    }

    ePrimeTry++;
  }

#ifdef HELIB_DEBUG
  std::cerr << "RecryptData::setAE(): e=" << e << ", e'=" << ePrime
            << std::endl;
#endif
}

bool RecryptData::operator==(const RecryptData& other) const
{
  if (mvec != other.mvec)
    return false;

  if (skHwt != other.skHwt)
    return false;

  return true;
}

// The main method
void RecryptData::init(const Context& context,
                       const NTL::Vec<long>& mvec_,
                       bool enableThick,
                       bool build_cache_,
                       bool minimal)
{
  if (alMod != nullptr) { // were we called for a second time?
    std::cerr << "@Warning: multiple calls to RecryptData::init\n";
    return;
  }

  // sanity check
  assertEq(computeProd(mvec_),
           context.getM(),
           "Cyclotomic polynomial mismatch");

  // Record the arguments to this function
  mvec = mvec_;
  build_cache = build_cache_;
  alsoThick = enableThick;

  bool mvec_ok = true;
  for (long i : range(mvec.length())) {
    NTL::Vec<NTL::Pair<long, long>> factors;
    factorize(factors, mvec[i]);
    if (factors.length() > 1)
      mvec_ok = false;
  }

  if (!mvec_ok) {
    Warning("prime power factorization recommended for bootstrapping");
  }

  skHwt = context.getHwt();
  e = context.getE();
  ePrime = context.getEPrime();

  long r = context.getAlMod().getR();

  // First part of Bootstrapping works wrt plaintext space p^{r'}
  alMod = std::make_shared<PAlgebraMod>(context.getZMStar(), e - ePrime + r);
  ea = std::make_shared<EncryptedArray>(context, *alMod);
  // Polynomial defaults to F0, PAlgebraMod explicitly given

  p2dConv = std::make_shared<PowerfulDCRT>(context, mvec);

  if (!enableThick)
    return;

  // Initialize the linear polynomial for unpacking the slots
  NTL::zz_pBak bak;
  bak.save();
  ea->getAlMod().restoreContext();
  long nslots = ea->size();
  long d = ea->getDegree();

  const NTL::Mat<NTL::zz_p>& CBi =
      ea->getDerived(PA_zz_p()).getNormalBasisMatrixInverse();

  std::vector<NTL::ZZX> LM;
  LM.resize(d);
  for (long i = 0; i < d; i++) // prepare the linear polynomial
    LM[i] = rep(CBi[i][0]);

  std::vector<NTL::ZZX> C;
  ea->buildLinPolyCoeffs(C, LM); // "build" the linear polynomial

  unpackSlotEncoding.resize(d); // encode the coefficients

  for (long j = 0; j < d; j++) {
    std::vector<NTL::ZZX> v(nslots);
    for (long k = 0; k < nslots; k++)
      v[k] = C[j];
    ea->encode(unpackSlotEncoding[j], v);
  }
  firstMap = std::make_shared<EvalMap>(*ea, minimal, mvec, true, build_cache);
  secondMap = std::make_shared<EvalMap>(context.getEA(),
                                        minimal,
                                        mvec,
                                        false,
                                        build_cache);
}

/********************************************************************/
/********************************************************************/

// Extract digits from fully packed slots
void extractDigitsPacked(Ctxt& ctxt,
                         long botHigh,
                         long r,
                         long ePrime,
                         const std::vector<NTL::ZZX>& unpackSlotEncoding, bool our_version = false, bool lazy = false);

// Extract digits from unpacked slots
void extractDigitsThin(Ctxt& ctxt, long botHigh, long r, long ePrime, bool our_version = false, bool lazy = false, std::vector<std::vector<long>> e_inner_compose_list = {{1}});

// Wrapper to make the above function part of the public libarary
void wrapExtractDigitsThin(Ctxt& ctxt, long botHigh, long r, bool our_version, bool lazy, std::vector<std::vector<long>> e_inner_compose_list, int nb_iterations) {
  Ctxt tmp(ctxt);
  long cap_before = ctxt.bitCapacity();
  int nb_relin_digit_extract;
  auto start = std::chrono::high_resolution_clock::now();
  auto duration_extract = start - start;

  /******************************
   ******** START LOOP **********
   ******************************/
  for (int iii = 0; iii < nb_iterations; iii++) {   // Average out multiple times
    ctxt = tmp;
    nb_relin_digit_extract = helib::nb_relin;
    start = std::chrono::high_resolution_clock::now();
    extractDigitsThin(ctxt, botHigh, r, /*ePrime*/ r, our_version, lazy, e_inner_compose_list);
    duration_extract += (std::chrono::high_resolution_clock::now() - start);
    nb_relin_digit_extract = helib::nb_relin - nb_relin_digit_extract;
  }
  /******************************
   ********* END LOOP ***********
   ******************************/

  // Make sure result is correct by taking negation (because homomorphic inner product is defined slightly differently in HElib)
  // Note that the correction for p equal to 2 is already done in the function extractDigitsThin
  ctxt.negate();

  //std::cout << "Number of relinearizations during digit extraction: " << nb_relin_digit_extract << std::endl;
  std::cout << "Number of digits to extract: " << botHigh << std::endl;
  std::cout << "Execution time" << std::endl;
  std::cout << "- Digit extraction: " << std::chrono::duration_cast<std::chrono::seconds>(duration_extract).count() / nb_iterations << " seconds." << std::endl;

  std::cout << "Noise capacity" << std::endl
            << "- Initial: " << cap_before << std::endl
            << "- Digit extract: " << cap_before - ctxt.bitCapacity() << std::endl;
}

// bootstrap a ciphertext to reduce noise
void PubKey::reCrypt(Ctxt& ctxt, bool our_version, bool lazy) const
{
  int nb_relin_digit_extract;
  long cap_first_map, cap_second_map, cap_digit_extract, cap_in_prod;
  auto start_time_bootstrapping = std::chrono::high_resolution_clock::now();
  auto total_time_digit_extract = start_time_bootstrapping - start_time_bootstrapping;

  HELIB_TIMER_START;

  // Some sanity checks for dummy ciphertext
  long ptxtSpace = ctxt.getPtxtSpace();
  if (ctxt.isEmpty())
    return;
  if (ctxt.parts.size() == 1 && ctxt.parts[0].skHandle.isOne()) {
    // Dummy encryption, just ensure that it is reduced mod p
    NTL::ZZX poly = to_ZZX(ctxt.parts[0]);
    for (long i = 0; i < poly.rep.length(); i++)
      poly[i] = NTL::to_ZZ(rem(poly[i], ptxtSpace));
    poly.normalize();
    ctxt.DummyEncrypt(poly);
    return;
  }

  // check that we have bootstrapping data
  assertTrue(recryptKeyID >= 0l, "No bootstrapping data");

  long p = getContext().getP();
  long r = getContext().getAlMod().getR();
  long p2r = getContext().getAlMod().getPPowR();

  long intFactor = ctxt.intFactor;

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  const RecryptData& rcData = getContext().getRcData();
  long e = rcData.e;
  long ePrime = rcData.ePrime;
  long p2ePrime = NTL::power_long(p, ePrime);
  long q = NTL::power_long(p, e) + 1;
  assertTrue(e >= r, "rcData.e must be at least alMod.r");

#ifdef HELIB_DEBUG
  std::cerr << "reCrypt: p=" << p << ", r=" << r << ", e=" << e
            << " ePrime=" << ePrime << ", q=" << q << std::endl;
  CheckCtxt(ctxt, "init");
#endif

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  assertEq(p2r % ptxtSpace, 0l, "ptxtSpace must divide p^r when bootstrapping");

  ctxt.dropSmallAndSpecialPrimes();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after mod down");
#endif

  HELIB_NTIMER_START(AAA_preProcess);

  // Make sure that this ciphertext is in canonical form
  if (!ctxt.inCanonicalForm())
    ctxt.reLinearize();

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / context.getSpecialPrimes();
  assertTrue(s <= context.getCtxtPrimes(), "prime set is messed up");
  if (s.card() > 3) { // leave only first three ciphertext primes
    long first = s.first();
    IndexSet s3(first, first + 2);
    s.retain(s3);
  }
  ctxt.modDownToSet(s);

  // key-switch to the bootstrapping key
  ctxt.reLinearize(recryptKeyID);

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after key switching");
#endif

  // "raw mod-switch" to the bootstrapping modulus q=p^e+1.
  std::vector<NTL::ZZX> zzParts; // the mod-switched parts, in ZZX format

  double mfac = ctxt.getContext().getZMStar().getNormBnd();
  double noise_est = ctxt.rawModSwitch(zzParts, q) * mfac;
  // noise_est is an upper bound on the L-infty norm of the scaled noise
  // in the pwrfl basis
  double noise_bnd =
      HELIB_MIN_CAP_FRAC * p2r * ctxt.getContext().boundForRecryption();
  // noise_bnd is the bound assumed in selecting the parameters
  double noise_rat = noise_est / noise_bnd;

  HELIB_STATS_UPDATE("raw-mod-switch-noise", noise_rat);

  if (noise_rat > 1) {
    // TODO: Turn the following preprocessor logics into a warnOrThrow function
    std::string message =
        "rawModSwitch scaled noise exceeds bound: " + std::to_string(noise_rat);
#ifdef HELIB_DEBUG
    Warning(message);
#else
    throw LogicError(message);
#endif
  }

  assertEq(zzParts.size(),
           (std::size_t)2,
           "Exactly 2 parts required for mod-switching in thin bootstrapping");

#ifdef HELIB_DEBUG
  if (dbgKey) {
    checkRecryptBounds(zzParts, dbgKey->getRecryptKey(), ctxt.getContext(), q);
  }
#endif

  std::vector<NTL::ZZX> v;
  v.resize(2);

  // Add multiples of q to make the zzParts divisible by p^{e'}
  for (long i : range(2)) {
    // make divisible by p^{e'}

    newMakeDivisible(zzParts[i], p2ePrime, q, ctxt.getContext(), v[i]);
  }

#ifdef HELIB_DEBUG
  if (dbgKey) {
    checkRecryptBounds_v(v, dbgKey->getRecryptKey(), ctxt.getContext(), q);
    checkCriticalValue(zzParts,
                       dbgKey->getRecryptKey(),
                       ctxt.getContext().getRcData(),
                       q);
  }
#endif

  for (long i : range(zzParts.size())) {
    zzParts[i] /= p2ePrime; // divide by p^{e'}
  }

  // NOTE: here we lose the intFactor associated with ctxt.
  // We will restore it below.
  ctxt = recryptEkey;

  ctxt.multByConstant(zzParts[1]);
  ctxt.addConstant(zzParts[0]);
  cap_in_prod = ctxt.bitCapacity();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after preProcess");
#endif
  HELIB_NTIMER_STOP(AAA_preProcess);

  // Move the powerful-basis coefficients to the plaintext slots
  HELIB_NTIMER_START(AAA_LinearTransform1);
  ctxt.getContext().getRcData().firstMap->apply(ctxt);
  HELIB_NTIMER_STOP(AAA_LinearTransform1);
  cap_first_map = ctxt.bitCapacity();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after LinearTransform1");
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  nb_relin_digit_extract = helib::nb_relin;
  auto start = std::chrono::high_resolution_clock::now();
  HELIB_NTIMER_START(AAA_extractDigitsPacked);
  extractDigitsPacked(ctxt,
                      e - ePrime,
                      r,
                      ePrime,
                      context.getRcData().unpackSlotEncoding, our_version, lazy);
  HELIB_NTIMER_STOP(AAA_extractDigitsPacked);
  total_time_digit_extract = std::chrono::high_resolution_clock::now() - start;
  nb_relin_digit_extract = helib::nb_relin - nb_relin_digit_extract;
  cap_digit_extract = ctxt.bitCapacity();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after extractDigitsPacked");
#endif

  // Move the slots back to powerful-basis coefficients
  HELIB_NTIMER_START(AAA_LinearTransform2);
  ctxt.getContext().getRcData().secondMap->apply(ctxt);
  HELIB_NTIMER_STOP(AAA_LinearTransform2);
  cap_second_map = ctxt.bitCapacity();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after linearTransform2");
#endif

  // restore intFactor
  if (intFactor != 1)
    ctxt.intFactor = NTL::MulMod(ctxt.intFactor, intFactor, ptxtSpace);

  //std::cout << "Number of relinearizations during digit extraction: " << nb_relin_digit_extract << std::endl;
  std::cout << "Number of digits to extract: " << e - ePrime << std::endl;
  std::cout << "Execution time" << std::endl;
  std::cout << "- Digit extraction: " << std::chrono::duration_cast<std::chrono::seconds>(total_time_digit_extract).count() << " seconds." << std::endl;
  std::cout << "- Total bootstrapping: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time_bootstrapping).count()
            << " seconds." << std::endl;

  std::cout << "Noise capacity" << std::endl
            << "- Initial: " << cap_in_prod << std::endl
            << "- Linear transformations: " << (cap_in_prod - cap_first_map) + (cap_digit_extract - cap_second_map) << std::endl
            << "- Digit extract: " << cap_first_map - cap_digit_extract << std::endl
            << "- Remaining: " << cap_second_map << std::endl;
}

#ifdef HELIB_BOOT_THREADS

// Extract digits from fully packed slots, multithreaded version
void extractDigitsPacked(Ctxt& ctxt,
                         long botHigh,
                         long r,
                         long ePrime,
                         const std::vector<NTL::ZZX>& unpackSlotEncoding, bool our_version, bool lazy)
{
  HELIB_TIMER_START;

  // Step 1: unpack the slots of ctxt
  HELIB_NTIMER_START(unpack);
  ctxt.cleanUp();

  // Apply the d automorphisms and store them in scratch area
  long d = ctxt.getContext().getOrdP();

  std::vector<Ctxt> unpacked(d, Ctxt(ZeroCtxtLike, ctxt));
  { // explicit scope to force all temporaries to be released
    std::vector<std::shared_ptr<DoubleCRT>> coeff_vector;
    std::vector<double> coeff_vector_sz;
    coeff_vector.resize(d);
    coeff_vector_sz.resize(d);

    HELIB_NTIMER_START(unpack1);
    for (long i = 0; i < d; i++) {
      coeff_vector[i] = std::make_shared<DoubleCRT>(unpackSlotEncoding[i],
                                                    ctxt.getContext(),
                                                    ctxt.getPrimeSet());
      coeff_vector_sz[i] = NTL::conv<double>(
          embeddingLargestCoeff(unpackSlotEncoding[i],
                                ctxt.getContext().getZMStar()));
    }
    HELIB_NTIMER_STOP(unpack1);

    HELIB_NTIMER_START(unpack2);
    std::vector<Ctxt> frob(d, Ctxt(ZeroCtxtLike, ctxt));

    NTL_EXEC_RANGE(d, first, last)
    // FIXME: implement using hoisting!
    for (long j = first; j < last; j++) { // process jth Frobenius
      frob[j] = ctxt;
      frob[j].frobeniusAutomorph(j);
      frob[j].cleanUp();
      // FIXME: not clear if we should call cleanUp here
    }
    NTL_EXEC_RANGE_END

    HELIB_NTIMER_STOP(unpack2);

    HELIB_NTIMER_START(unpack3);
    Ctxt tmp1(ZeroCtxtLike, ctxt);
    for (long i = 0; i < d; i++) {
      for (long j = 0; j < d; j++) {
        tmp1 = frob[j];
        tmp1.multByConstant(*coeff_vector[mcMod(i + j, d)],
                            coeff_vector_sz[mcMod(i + j, d)]);
        unpacked[i] += tmp1;
      }
    }
    HELIB_NTIMER_STOP(unpack3);
  }
  HELIB_NTIMER_STOP(unpack);

  //#ifdef HELIB_DEBUG
  //  CheckCtxt(unpacked[0], "after unpack");
  //#endif

  NTL_EXEC_RANGE(d, first, last)
  for (long i = first; i < last; i++) {
    extractDigitsThin(unpacked[i], botHigh, r, ePrime, our_version, lazy);
  }
  NTL_EXEC_RANGE_END

  //#ifdef HELIB_DEBUG
  // CheckCtxt(unpacked[0], "before repack");
  //#endif

  // Step 3: re-pack the slots
  HELIB_NTIMER_START(repack);
  const EncryptedArray& ea2 = ctxt.getContext().getEA();
  NTL::ZZX xInSlots;
  std::vector<NTL::ZZX> xVec(ea2.size());
  ctxt = unpacked[0];
  for (long i = 1; i < d; i++) {
    x2iInSlots(xInSlots, i, xVec, ea2);
    unpacked[i].multByConstant(xInSlots);
    ctxt += unpacked[i];
  }
  HELIB_NTIMER_STOP(repack);
  //#ifdef HELIB_DEBUG
  // CheckCtxt(ctxt, "after repack");
  //#endif
}

#else

// Extract digits from fully packed slots
void extractDigitsPacked(Ctxt& ctxt,
                         long botHigh,
                         long r,
                         long ePrime,
                         const std::vector<NTL::ZZX>& unpackSlotEncoding, bool our_version, bool lazy)
{
  HELIB_TIMER_START;

  // Step 1: unpack the slots of ctxt
  HELIB_NTIMER_START(unpack);
  ctxt.cleanUp();

  // Apply the d automorphisms and store them in scratch area
  long d = ctxt.getContext().getOrdP();

  std::vector<Ctxt> unpacked(d, Ctxt(ZeroCtxtLike, ctxt));
  { // explicit scope to force all temporaries to be released
    std::vector<std::shared_ptr<DoubleCRT>> coeff_vector;
    std::vector<double> coeff_vector_sz;
    coeff_vector.resize(d);
    coeff_vector_sz.resize(d);
    for (long i = 0; i < d; i++) {
      coeff_vector[i] = std::make_shared<DoubleCRT>(unpackSlotEncoding[i],
                                                    ctxt.getContext(),
                                                    ctxt.getPrimeSet());
      coeff_vector_sz[i] = NTL::conv<double>(
          embeddingLargestCoeff(unpackSlotEncoding[i],
                                ctxt.getContext().getZMStar()));
    }

    Ctxt tmp1(ZeroCtxtLike, ctxt);
    Ctxt tmp2(ZeroCtxtLike, ctxt);

    // FIXME: implement using hoisting!
    for (long j = 0; j < d; j++) { // process jth Frobenius
      tmp1 = ctxt;
      tmp1.frobeniusAutomorph(j);
      tmp1.cleanUp();
      // FIXME: not clear if we should call cleanUp here

      for (long i = 0; i < d; i++) {
        tmp2 = tmp1;
        tmp2.multByConstant(*coeff_vector[mcMod(i + j, d)],
                            coeff_vector_sz[mcMod(i + j, d)]);
        unpacked[i] += tmp2;
      }
    }
  }
  HELIB_NTIMER_STOP(unpack);

  //#ifdef HELIB_DEBUG
  //  CheckCtxt(unpacked[0], "after unpack");
  //#endif

  for (long i = 0; i < (long)unpacked.size(); i++) {
    extractDigitsThin(unpacked[i], botHigh, r, ePrime, our_version, lazy);
  }

  //#ifdef HELIB_DEBUG
  //  CheckCtxt(unpacked[0], "before repack");
  //#endif

  // Step 3: re-pack the slots
  HELIB_NTIMER_START(repack);
  const EncryptedArray& ea2 = ctxt.getContext().getEA();
  NTL::ZZX xInSlots;
  std::vector<NTL::ZZX> xVec(ea2.size());
  ctxt = unpacked[0];
  for (long i = 1; i < d; i++) {
    x2iInSlots(xInSlots, i, xVec, ea2);
    unpacked[i].multByConstant(xInSlots);
    ctxt += unpacked[i];
  }
  HELIB_NTIMER_STOP(repack);
}

#endif

// Use packed bootstrapping, so we can bootstrap all in just one go.
void packedRecrypt(const CtPtrs& cPtrs,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea)
{
  PubKey& pKey = (PubKey&)cPtrs[0]->getPubKey();

  // Allocate temporary ciphertexts for the recryption
  int nPacked = divc(cPtrs.size(), ea.getDegree()); // ceil(totalNum/d)
  std::vector<Ctxt> cts(nPacked, Ctxt(pKey));

  repack(CtPtrs_vectorCt(cts), cPtrs, ea); // pack ciphertexts
  //  cout << "@"<< lsize(cts)<<std::flush;
  for (Ctxt& c : cts) {   // then recrypt them
    c.reducePtxtSpace(2); // we only have recryption data for binary ctxt
    pKey.reCrypt(c);
  }
  unpack(cPtrs, CtPtrs_vectorCt(cts), ea, unpackConsts);
}

// recrypt all ctxt at level < belowLvl
void packedRecrypt(const CtPtrs& array,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea,
                   long belowLvl)
{
  std::vector<Ctxt*> v;
  for (long i = 0; i < array.size(); i++)
    if (array.isSet(i) && !array[i]->isEmpty() &&
        array[i]->bitCapacity() < belowLvl * (array[i]->getContext().BPL()))
      v.push_back(array[i]);
  packedRecrypt(CtPtrs_vectorPt(v), unpackConsts, ea);
}
void packedRecrypt(const CtPtrMat& m,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea,
                   long belowLvl)
{
  std::vector<Ctxt*> v;
  for (long i = 0; i < m.size(); i++)
    for (long j = 0; j < m[i].size(); j++)
      if (m[i].isSet(j) && !m[i][j]->isEmpty() &&
          m[i][j]->bitCapacity() < belowLvl * (m[i][j]->getContext().BPL()))
        v.push_back(m[i][j]);
  packedRecrypt(CtPtrs_vectorPt(v), unpackConsts, ea);
}

//===================== Thin Bootstrapping stuff ==================

void ThinRecryptData::init(const Context& context,
                           const NTL::Vec<long>& mvec_,
                           bool alsoThick,
                           bool build_cache_,
                           bool minimal)
{
  RecryptData::init(context, mvec_, alsoThick, build_cache_, minimal);
  coeffToSlot =
      std::make_shared<ThinEvalMap>(*ea, minimal, mvec, true, build_cache);
  slotToCoeff = std::make_shared<ThinEvalMap>(context.getEA(),
                                              minimal,
                                              mvec,
                                              false,
                                              build_cache);
}

// Extract digits from thinly packed slots

long fhe_force_chen_han = 0;

// Evaluate the optimized bit extraction polynomials using the multivariate strategy (only for p = 2 and e <= 16)
// We always use the same set of polynomials, regardless of e and how many results are actually required
// Parameters:
// - Ciphertext in which to evaluate the polynomials
// - Result vector that will include ciphertexts and precisions
// - Size of the row: distance between leftmost and rightmost digit (counting is done from input to output)
void rowComputationMultivariate(const Ctxt& ctxt, std::vector<std::pair<Ctxt, long>>& ctxtEval, long rowSize) {
    // List of polynomials:
    // * f2 = x^2
    // * f4 = (f2)^2
    // * f8 = 112*f2 + (94*f2 + 121*f4)^2
    // * f16 = 11136*f4 - (15364*f4 - 14115*f8) * (28504*f2 + 8968*f4 - f8)

    // Evaluate the above polynomials
    Ctxt tmp1(ZeroCtxtLike, ctxt);
    Ctxt tmp2(ZeroCtxtLike, ctxt);
    Ctxt tmp3(ZeroCtxtLike, ctxt);
    long init_length = ctxtEval.size();
    if (rowSize >= 2) {
        ctxtEval.push_back(std::pair<Ctxt, long>(ctxt, 2));
        std::get<0>(ctxtEval[init_length + 0]).square();
    }
    if (rowSize >= 3) {
        ctxtEval.push_back(std::pair<Ctxt, long>(std::get<0>(ctxtEval[init_length + 0]), 4));
        std::get<0>(ctxtEval[init_length + 1]).square();
    }
    if (rowSize >= 5) {
        ctxtEval.push_back(std::pair<Ctxt, long>(std::get<0>(ctxtEval[init_length + 0]), 8));
        std::get<0>(ctxtEval[init_length + 2]).multByConstant((long) 112);

        tmp1 = std::get<0>(ctxtEval[init_length + 0]);
        tmp1.multByConstant((long) 94);

        tmp2 = std::get<0>(ctxtEval[init_length + 1]);
        tmp2.multByConstant((long) 121);
        tmp1.addCtxt(tmp2);

        tmp1.square();

        std::get<0>(ctxtEval[init_length + 2]).addCtxt(tmp1);
    }
    if (rowSize >= 9) {
        ctxtEval.push_back(std::pair<Ctxt, long>(std::get<0>(ctxtEval[init_length + 1]), 16));
        std::get<0>(ctxtEval[init_length + 3]).multByConstant((long) 11136);

        tmp1 = std::get<0>(ctxtEval[init_length + 1]);
        tmp1.multByConstant((long) 15364);

        tmp2 = std::get<0>(ctxtEval[init_length + 2]);
        tmp2.multByConstant((long) 14115);
        tmp1.addCtxt(tmp2, true);

        tmp2 = std::get<0>(ctxtEval[init_length + 0]);
        tmp2.multByConstant((long) 28504);

        tmp3 = std::get<0>(ctxtEval[init_length + 1]);
        tmp3.multByConstant((long) 8968);
        tmp2.addCtxt(tmp3);

        tmp2.addCtxt(std::get<0>(ctxtEval[init_length + 2]), true);

        tmp1.multiplyBy(tmp2);

        std::get<0>(ctxtEval[init_length + 3]).addCtxt(tmp1, true);
    }
}

// Evaluate the optimized digit extraction polynomials using the even/odd strategy (this optimization is not
// done when p = 2 and e_inner > 1) and the function composition approach
// Parameters:
// - Ciphertext in which to evaluate the polynomials
// - Result vector that will include ciphertexts and precisions
// - Size of the inner part of the trapezoid (the v x v triangle): only the current row is considered
// - Size of the row: distance between leftmost and rightmost digit (counting is done from input to output)
// - Flag for lazy relinearization or default strategy
// - Precision exponent of the input ciphertext (relevant for function composition approach)
void rowComputationComposition(const Ctxt& ctxt, std::vector<std::pair<Ctxt, long>>& ctxtEval, long triangleSize, long rowSize, bool lazy, long e_inner) {
    // Find index of primes in list
    long index = nb_primes;
    for (long iii = 0; iii < nb_primes; iii++)
        if (primes_list[iii] == ctxt.getContext().getP())
            index = iii;
    assertTrue(index < nb_primes, "No polynomials generated for given prime.");

    // Find index of e_inner in list
    long inner_index = nb_e_inner;
    for (long iii = 0; iii < nb_e_inner; iii++)
        if (e_inner_list[iii] == e_inner)
            inner_index = iii;
    assertTrue(inner_index < nb_e_inner, "No polynomials generated for given e_inner.");

    // Always minimize the multiplicative depth as a rule of thumb (see paper for description of heuristic)
    std::vector<NTL::ZZX> polynomials;  // Polynomials to evaluate
    std::vector<long> precisions;       // Precisions of the polynomials
    for (long exponent = 1; e_inner * pow(2, exponent - 1) < triangleSize; exponent++) {  // As long as previous iteration did not reach triangle size
        long precision = e_inner * pow(2, exponent);
        if (precision > triangleSize) { // If we are outside small triangle already
            if (precision < rowSize)    // If we did not yet reach the entire row size
                precision = triangleSize;
            else
                precision = rowSize;
        }
        assertTrue(precision - e_inner - 1 < (long)polynomial_vector[index][inner_index].size(), "Not sufficiently many polynomials generated for the given prime.");
        polynomials.push_back(polynomial_vector[index][inner_index][precision - e_inner - 1]);  // Polynomials are loaded from precision = e_inner + 1
        precisions.push_back(precision);
    }

    // Possibly add one more polynomial
    if (precisions.empty() || precisions.back() < rowSize) {
        assertTrue(rowSize - e_inner - 1 < (long)polynomial_vector[index][inner_index].size(), "Not sufficiently many polynomials generated for the given prime.");
        polynomials.push_back(polynomial_vector[index][inner_index][rowSize - e_inner - 1]);
        precisions.push_back(rowSize);
    }

    // Evaluate polynomials using Paterson-Stockmeyer
    std::vector<Ctxt> result;
    customPolyEval(result, polynomials, ctxt, lazy);

    // Put result in return argument
    for (long index = 0; index < (long)result.size(); index++) {
        ctxtEval.push_back(std::pair<Ctxt, long>(result[index], precisions[index]));
    }
}

// Evaluate the optimized digit extraction polynomials
// Same functionality as two functions above, except that we pass a list of values for e_inner: this gives the different splitting
// values, where the first one indicates the precision of the input ciphertext (should normally be 1)
void rowComputationGeneral(const Ctxt& ctxt, std::vector<std::pair<Ctxt, long>>& ctxtEval, long triangleSize, long rowSize, bool lazy, std::vector<long> e_inner_compose_list) {
    e_inner_compose_list.push_back(rowSize);
    ctxtEval.push_back(std::pair<Ctxt, long>(ctxt, e_inner_compose_list.front()));

    // Call the usual function for evaluating digit extraction polynomials
    for (long index = 1; index < (long)e_inner_compose_list.size(); index++) {
        long e_inner_previous = e_inner_compose_list[index - 1];
        long e_inner = e_inner_compose_list[index];
        if ((ctxt.getContext().getP() == 2) && (e_inner_previous == 1) && (e_inner <= 16))
            rowComputationMultivariate(std::get<0>(ctxtEval.back()), ctxtEval, std::min(rowSize, e_inner));
        else
            rowComputationComposition(std::get<0>(ctxtEval.back()), ctxtEval, std::min(triangleSize, e_inner), std::min(rowSize, e_inner), lazy, e_inner_previous);
    }
}

// Our improved digit extraction algorithm
void customExtractDigitsThin(Ctxt& ctxt, long botHigh, long r, bool lazy, std::vector<std::vector<long>> e_inner_compose_list) {
    // Apply correction for p = 2, because balanced digit representation does not exist
    if (ctxt.getContext().getP() == 2)
        ctxt.addConstant(lround(pow(ctxt.getContext().getP(), botHigh) / 2));

    // Format of the vectors below: keep pair of ciphertext and number, where
    // the number indicates the exponent to which the ciphertext is defined
    // - For number e, the corresponding ciphertext is defined mod p^e (so
    //   the e lower digits are correct and the other ones are garbage)
    // - The numbers are guaranteed to be in increasing order
    std::vector<std::pair<Ctxt, long>> ctxtRows(botHigh, std::make_pair(ctxt, botHigh + r));  // Represent each row of the trapezoid (not including result)
    for (int row = 0; row < botHigh; row++) {
        // Evaluate necessary polynomials only
        std::vector<std::pair<Ctxt, long>> ctxtEval;                                          // Store evaluation of digit extraction polynomials
        rowComputationGeneral(std::get<0>(ctxtRows[row]), ctxtEval, botHigh - row, botHigh + r - row, lazy, e_inner_compose_list[std::min(row, (int)e_inner_compose_list.size() - 1)]);

        // Determine starting values for next rows based on the necessary precision
        for (int nextRow = row + 1; nextRow < botHigh; nextRow++) {   // Update next rows with the result from above
            // Check if we already have result with required precision (not possible for row + 1)
            if ((nextRow > row + 1) && (std::get<1>(ctxtRows[nextRow - 1]) + row + 1 >= nextRow + 1)) { // Compare precisions (interpret them wrt highest exponent botHigh + r)
                ctxtRows[nextRow] = ctxtRows[nextRow - 1];
            } else {
                // Loop over the result from polynomial evaluation
                for (std::pair<Ctxt, long> tup : ctxtEval) {
                    if (std::get<1>(tup) + row >= nextRow + 1) {                        // Compare precisions (interpret them wrt highest exponent botHigh + r)
                        std::get<0>(ctxtRows[nextRow]).addCtxt(std::get<0>(tup), true); // Subtract extracted digit
                        std::get<0>(ctxtRows[nextRow]).divideByP();                     // Divide by p
                        std::get<1>(ctxtRows[nextRow]) = std::min(std::get<1>(ctxtRows[nextRow]), std::get<1>(tup)) - 1;  // Update stored precision
                        break;
                    }
                }
            }
        }

        // Finally compute the result in a similar way as above
        // Check if we already have result with required precision (not possible for last row)
        if ((botHigh > row + 1) && (std::get<1>(ctxtRows.back()) + row + 1 >= botHigh + r)) { // Compare precisions (interpret them wrt highest exponent botHigh + r)
            ctxt = std::get<0>(ctxtRows.back());
        } else {
            ctxt.addCtxt(std::get<0>(ctxtEval.back()), true); // Subtract extracted digit
            ctxt.divideByP();                                 // Divide by p
        }
    }

    // Necessary due to different version of homomorphic inner product in HElib
    ctxt.negate();
}

// Built-in digit extraction algorithm (we just call our own function inside)
void extractDigitsThin(Ctxt& ctxt, long botHigh, long r, long ePrime, bool our_version, bool lazy, std::vector<std::vector<long>> e_inner_compose_list)
{
  HELIB_TIMER_START;

  if (ePrime < r) {     // For HElib version of homomorphic inner product, this is not allowed if we want to extract upper digits only
    std::cout << "Warning: unfortunate choice of parameters (complexity of digit extraction is unnecessarily high because e' < r).\nPatermeters: e' = "
              << ePrime << " and r = " << r << std::endl;
    if (our_version)
        throw RuntimeError("Bad parameter choice. See warning above.");
  }

  // Call our own digit extraction function
  if (our_version) {
    customExtractDigitsThin(ctxt, botHigh, r, lazy, e_inner_compose_list);
  } else {
    std::cout << "Starting built-in digit extraction." << std::endl;
    Ctxt unpacked(ctxt);
    unpacked.cleanUp();

    std::vector<Ctxt> scratch;

    long p = ctxt.getContext().getP();
    long p2r = NTL::power_long(p, r);
    long topHigh = botHigh + r - 1;

    // degree Chen/Han technique is p^{bot-1}(p-1)r
    // degree of basic technique is p^{bot-1}p^r,
    //     or p^{bot-1}p^{r-1} if p==2, r > 1, and bot+r > 2

    bool use_chen_han = false;
    if (r > 1) {
      double chen_han_cost = log(p - 1) + log(r);
      double basic_cost;
      if (p == 2 && r > 2 && botHigh + r > 2)
        basic_cost = (r - 1) * log(p);
      else
        basic_cost = r * log(p);

      // std::cerr << "*** basic: " << basic_cost << "\n";
      // std::cerr << "*** chen/han: " << chen_han_cost << "\n";

      double thresh = 1.5;
      if (p == 2)
        thresh = 1.75;
      // increasing thresh makes chen_han less likely to be chosen.
      // For p == 2, the basic algorithm is just squaring,
      // and so is a bit cheaper, so we raise thresh a bit.
      // This is all a bit heuristic.

      if (basic_cost > thresh * chen_han_cost)
        use_chen_han = true;
    }

    if (fhe_force_chen_han > 0)
      use_chen_han = true;
    else if (fhe_force_chen_han < 0)
      use_chen_han = false;

    if (!use_chen_han)
      std::cout << "Warning: Chen/Han digit extraction is not being used." << std::endl;

    if (use_chen_han) {
      // use Chen and Han technique

      extendExtractDigits(scratch, unpacked, botHigh, r);

#if 0
      for (long i: range(scratch.size())) {
        CheckCtxt(scratch[i], "**");
      }
#endif

      for (long j = 0; j < botHigh; j++) {
        unpacked -= scratch[j];
        unpacked.divideByP();
      }

      if (p == 2 && botHigh > 0) // For p==2, subtract also the previous bit
        unpacked += scratch[botHigh - 1];
      unpacked.negate();

      if (r > ePrime) { // Add in digits from the bottom part, if any
        long topLow = r - 1 - ePrime;
        Ctxt tmp = scratch[topLow];
        for (long j = topLow - 1; j >= 0; --j) {
          tmp.multByP();
          tmp += scratch[j];
        }
        if (ePrime > 0)
          tmp.multByP(ePrime); // multiply by p^e'
        unpacked += tmp;
      }
      unpacked.reducePtxtSpace(p2r); // Our plaintext space is now mod p^r

      ctxt = unpacked;
    } else {

      if (p == 2 && r > 2 && topHigh + 1 > 2)
        topHigh--; // For p==2 we sometime get a bit for free

      extractDigits(scratch, unpacked, topHigh + 1);

      // set unpacked = -\sum_{j=botHigh}^{topHigh} scratch[j] * p^{j-botHigh}
      if (topHigh >= LONG(scratch.size())) {
        topHigh = scratch.size() - 1;
        std::cerr << " @ suspect: not enough digits in extractDigitsPacked\n";
      }

      unpacked = scratch[topHigh];
      for (long j = topHigh - 1; j >= botHigh; --j) {
        unpacked.multByP();
        unpacked += scratch[j];
      }
      if (p == 2 && botHigh > 0) // For p==2, subtract also the previous bit
        unpacked += scratch[botHigh - 1];
      unpacked.negate();

      if (r > ePrime) { // Add in digits from the bottom part, if any
        long topLow = r - 1 - ePrime;
        Ctxt tmp = scratch[topLow];
        for (long j = topLow - 1; j >= 0; --j) {
          tmp.multByP();
          tmp += scratch[j];
        }
        if (ePrime > 0)
          tmp.multByP(ePrime); // multiply by p^e'
        unpacked += tmp;
      }
      unpacked.reducePtxtSpace(p2r); // Our plaintext space is now mod p^r
      ctxt = unpacked;
    }
  }
}

// Hack to get at private fields of public key
struct PubKeyHack
{                         // The public key
  const Context& context; // The context

  //! @var Ctxt pubEncrKey
  //! The public encryption key is an encryption of 0,
  //! relative to the first secret key
  Ctxt pubEncrKey;

  std::vector<long> skHwts;            // The Hamming weight of the secret keys
  std::vector<KeySwitch> keySwitching; // The key-switching matrices

  // The keySwitchMap structure contains pointers to key-switching matrices
  // for re-linearizing automorphisms. The entry keySwitchMap[i][n] contains
  // the index j such that keySwitching[j] is the first matrix one needs to
  // use when re-linearizing s_i(X^n).
  std::vector<std::vector<long>> keySwitchMap;

  NTL::Vec<int> KS_strategy; // NTL Vec's support I/O, which is
                             // more convenient

  // bootstrapping data

  long recryptKeyID; // index of the bootstrapping key
  Ctxt recryptEkey;  // the key itself, encrypted under key #0
};

// bootstrap a ciphertext to reduce noise
void PubKey::thinReCrypt(Ctxt& ctxt, bool our_version, bool lazy, int nb_iterations) const
{
  int nb_relin_digit_extract = 0;
  long cap_first_map = 0, cap_second_map = 0, cap_digit_extract = 0, cap_in_prod = 0, cap_start = 0;
  auto start_time_bootstrapping = std::chrono::high_resolution_clock::now();
  auto total_time_digit_extract = start_time_bootstrapping - start_time_bootstrapping;
  cap_start = ctxt.bitCapacity();

  // Declare variables outside loop already so that we can print them at the end
  long e = 0, ePrime = 0;

  /******************************
   ******** START LOOP **********
   ******************************/
  for (int iii = 0; iii < nb_iterations; iii++) {  // Average out multiple times
  HELIB_TIMER_START;

  // Some sanity checks for dummy ciphertext
  long ptxtSpace = ctxt.getPtxtSpace();
  if (ctxt.isEmpty())
    return;

  if (ctxt.parts.size() == 1 && ctxt.parts[0].skHandle.isOne()) {
    // Dummy encryption, just ensure that it is reduced mod p
    NTL::ZZX poly = to_ZZX(ctxt.parts[0]);
    for (long i = 0; i < poly.rep.length(); i++)
      poly[i] = NTL::to_ZZ(rem(poly[i], ptxtSpace));
    poly.normalize();
    ctxt.DummyEncrypt(poly);
    return;
  }

  // check that we have bootstrapping data
  assertTrue(recryptKeyID >= 0l, "Bootstrapping data not present");

  long p = ctxt.getContext().getP();
  long r = ctxt.getContext().getAlMod().getR();
  long p2r = ctxt.getContext().getAlMod().getPPowR();

  long intFactor = ctxt.intFactor;

  const ThinRecryptData& trcData = ctxt.getContext().getRcData();

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  e = trcData.e;
  ePrime = trcData.ePrime;
  long p2ePrime = NTL::power_long(p, ePrime);
  long q = NTL::power_long(p, e) + 1;
  assertTrue(e >= r, "trcData.e must be at least alMod.r");

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  assertEq(p2r % ptxtSpace,
           0l,
           "ptxtSpace must divide p^r when thin bootstrapping");

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "init");
#endif

  ctxt.dropSmallAndSpecialPrimes();

#define DROP_BEFORE_THIN_RECRYPT
#define THIN_RECRYPT_NLEVELS (3)
#ifdef DROP_BEFORE_THIN_RECRYPT
  // experimental code...we should drop down to a reasonably low level
  // before doing the first linear map.
  long first = context.getCtxtPrimes().first();
  long last = std::min(context.getCtxtPrimes().last(),
                       first + THIN_RECRYPT_NLEVELS - 1);
  ctxt.bringToSet(IndexSet(first, last));
#endif

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after mod down");
#endif

  // Move the slots to powerful-basis coefficients
  HELIB_NTIMER_START(AAA_slotToCoeff);
  trcData.slotToCoeff->apply(ctxt);
  HELIB_NTIMER_STOP(AAA_slotToCoeff);
  cap_first_map = ctxt.bitCapacity();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after slotToCoeff");
#endif

  HELIB_NTIMER_START(AAA_bootKeySwitch);

  // Make sure that this ciphertext is in canonical form
  if (!ctxt.inCanonicalForm())
    ctxt.reLinearize();

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / context.getSpecialPrimes();
  assertTrue(s <= context.getCtxtPrimes(), "prime set is messed up");
  if (s.card() > 3) { // leave only first three ciphertext primes
    long first = s.first();
    IndexSet s3(first, first + 2);
    s.retain(s3);
  }
  ctxt.modDownToSet(s);

  // key-switch to the bootstrapping key
  ctxt.reLinearize(recryptKeyID);

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after key switching");
#endif

  // "raw mod-switch" to the bootstrapping mosulus q=p^e+1.
  std::vector<NTL::ZZX> zzParts; // the mod-switched parts, in ZZX format

  double mfac = ctxt.getContext().getZMStar().getNormBnd();
  double noise_est = ctxt.rawModSwitch(zzParts, q) * mfac;
  // noise_est is an upper bound on the L-infty norm of the scaled noise
  // in the pwrfl basis
  double noise_bnd =
      HELIB_MIN_CAP_FRAC * p2r * ctxt.getContext().boundForRecryption();
  // noise_bnd is the bound assumed in selecting the parameters
  double noise_rat = noise_est / noise_bnd;

  HELIB_STATS_UPDATE("raw-mod-switch-noise", noise_rat);

  if (noise_rat > 1) {
    // TODO: Turn the following preprocessor logics into a warnOrThrow function
    std::string message =
        "rawModSwitch scaled noise exceeds bound: " + std::to_string(noise_rat);
#ifdef HELIB_DEBUG
    Warning(message);
#else
    throw LogicError(message);
#endif
  }

  assertEq(zzParts.size(),
           (std::size_t)2,
           "Exactly 2 parts required for mod-switching in thin bootstrapping");

#ifdef HELIB_DEBUG
  if (dbgKey) {
    checkRecryptBounds(zzParts, dbgKey->getRecryptKey(), ctxt.getContext(), q);
  }
#endif

  std::vector<NTL::ZZX> v;
  v.resize(2);

  // Add multiples of q to make the zzParts divisible by p^{e'}
  for (long i : range(2)) {
    // make divisible by p^{e'}

    newMakeDivisible(zzParts[i], p2ePrime, q, ctxt.getContext(), v[i]);
  }

#ifdef HELIB_DEBUG
  if (dbgKey) {
    checkRecryptBounds_v(v, dbgKey->getRecryptKey(), ctxt.getContext(), q);
    checkCriticalValue(zzParts,
                       dbgKey->getRecryptKey(),
                       ctxt.getContext().getRcData(),
                       q);
  }
#endif

  for (long i : range(zzParts.size())) {
    zzParts[i] /= p2ePrime; // divide by p^{e'}
  }

  // NOTE: here we lose the intFactor associated with ctxt.
  // We will restore it below.
  ctxt = recryptEkey;

  ctxt.multByConstant(zzParts[1]);
  ctxt.addConstant(zzParts[0]);
  cap_in_prod = ctxt.bitCapacity();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after bootKeySwitch");
#endif

  HELIB_NTIMER_STOP(AAA_bootKeySwitch);

  // Move the powerful-basis coefficients to the plaintext slots
  HELIB_NTIMER_START(AAA_coeffToSlot);
  trcData.coeffToSlot->apply(ctxt);
  HELIB_NTIMER_STOP(AAA_coeffToSlot);
  cap_second_map = ctxt.bitCapacity();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after coeffToSlot");
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  nb_relin_digit_extract = helib::nb_relin;
  auto start = std::chrono::high_resolution_clock::now();
  HELIB_NTIMER_START(AAA_extractDigitsThin);
  extractDigitsThin(ctxt, e - ePrime, r, ePrime, our_version, lazy);
  HELIB_NTIMER_STOP(AAA_extractDigitsThin);
  total_time_digit_extract += (std::chrono::high_resolution_clock::now() - start);
  nb_relin_digit_extract = helib::nb_relin - nb_relin_digit_extract;
  cap_digit_extract = ctxt.bitCapacity();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after extractDigitsThin");
#endif

  // restore intFactor
  if (intFactor != 1)
    ctxt.intFactor = NTL::MulMod(ctxt.intFactor, intFactor, ptxtSpace);
  }
  /******************************
   ********* END LOOP ***********
   ******************************/

  //std::cout << "Number of relinearizations during digit extraction: " << nb_relin_digit_extract << std::endl;
  std::cout << "Number of digits to extract: " << e - ePrime << std::endl;
  std::cout << "Execution time" << std::endl;
  std::cout << "- Digit extraction: " << std::chrono::duration_cast<std::chrono::seconds>(total_time_digit_extract).count() / nb_iterations << " seconds." << std::endl;
  std::cout << "- Total bootstrapping: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time_bootstrapping).count() / nb_iterations
            << " seconds." << std::endl;

  std::cout << "Noise capacity" << std::endl
            << "- Initial: " << cap_in_prod << std::endl
            << "- Linear transformations: " << (cap_start - cap_first_map) + (cap_in_prod - cap_second_map) << std::endl
            << "- Digit extract: " << cap_second_map - cap_digit_extract << std::endl
            << "- Remaining: " << cap_digit_extract - (cap_start - cap_first_map) << std::endl;
}

#ifdef HELIB_DEBUG

static void checkCriticalValue(const std::vector<NTL::ZZX>& zzParts,
                               const DoubleCRT& sKey,
                               const RecryptData& rcData,
                               long q)
{
  NTL::ZZX ptxt;
  rawDecrypt(ptxt, zzParts, sKey); // no mod q

  NTL::Vec<NTL::ZZ> powerful;
  rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
  NTL::xdouble max_pwrfl = NTL::conv<NTL::xdouble>(largestCoeff(powerful));
  double critical_value = NTL::conv<double>((max_pwrfl / q) / q);

  vecRed(powerful, powerful, q, false);
  max_pwrfl = NTL::conv<NTL::xdouble>(largestCoeff(powerful));
  critical_value += NTL::conv<double>(max_pwrfl / q);

  HELIB_STATS_UPDATE("critical-value", critical_value);

  std::cerr << "=== critical_value=" << critical_value;
  if (critical_value > 0.5)
    std::cerr << " BAD-BOUND";

  std::cerr << "\n";
}

static void checkRecryptBounds(const std::vector<NTL::ZZX>& zzParts,
                               const DoubleCRT& sKey,
                               const Context& context,
                               long q)
{
  const RecryptData& rcData = context.getRcData();
  double coeff_bound = context.boundForRecryption();
  long p2r = context.getAlMod().getPPowR();

  NTL::ZZX ptxt;
  rawDecrypt(ptxt, zzParts, sKey); // no mod q

  NTL::Vec<NTL::ZZ> powerful;
  rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
  double max_pwrfl = NTL::conv<double>(largestCoeff(powerful));
  double ratio = max_pwrfl / (2 * q * coeff_bound);

  HELIB_STATS_UPDATE("|x|/bound", ratio);

  std::cerr << "=== |x|/bound=" << ratio;
  if (ratio > 1.0)
    std::cerr << " BAD-BOUND";

  vecRed(powerful, powerful, q, false);
  max_pwrfl = NTL::conv<double>(largestCoeff(powerful));
  ratio = max_pwrfl / (2 * p2r * coeff_bound);

  HELIB_STATS_UPDATE("|x%q|/bound", ratio);

  std::cerr << ", (|x%q|)/bound=" << ratio;
  if (ratio > 1.0)
    std::cerr << " BAD-BOUND";

  std::cerr << "\n";
}

static void checkRecryptBounds_v(const std::vector<NTL::ZZX>& v,
                                 const DoubleCRT& sKey,
                                 const Context& context,
                                 UNUSED long q)
{
  const RecryptData& rcData = context.getRcData();

  long p = context.getP();
  long e = rcData.e;
  long p2e = NTL::power_long(p, e);
  long ePrime = rcData.ePrime;
  long p2ePrime = NTL::power_long(p, ePrime);
  long phim = context.getPhiM();

  double fudge = compute_fudge(p2ePrime, p2e);

  double coeff_bound = context.boundForRecryption() * fudge;

  double sigma = context.stdDevForRecryption() * fudge;

  NTL::ZZX ptxt;
  rawDecrypt(ptxt, v, sKey); // no mod q

  NTL::Vec<NTL::ZZ> powerful;
  rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
  double max_pwrfl = NTL::conv<double>(largestCoeff(powerful));

  double denom = p2ePrime * coeff_bound;
  double ratio = max_pwrfl / denom;

  HELIB_STATS_UPDATE("|v|/bound", ratio);

  std::cerr << "=== |v|/bound=" << ratio;
  if (ratio > 1.0)
    std::cerr << " BAD-BOUND";
  std::cerr << "\n";

  ptxt -= v[0]; // so now ptxt is just sKey * v[1]
  rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);

  assertEq(powerful.length(), phim, "length should be phim");

  double ran_pwrfl = NTL::conv<double>(powerful[NTL::RandomBnd(phim)]);
  // pick a random coefficient in the poweful basis

  double std_devs = fabs(ran_pwrfl) / (p2ePrime * sigma);
  // number of standard deviations away from mean

  // update various indicator variables
  HELIB_STATS_UPDATE("sigma_0_5", double(std_devs <= 0.5)); // 0.383
  HELIB_STATS_UPDATE("sigma_1_0", double(std_devs <= 1.0)); // 0.683
  HELIB_STATS_UPDATE("sigma_1_5", double(std_devs <= 1.5)); // 0.866
  HELIB_STATS_UPDATE("sigma_2_0", double(std_devs <= 2.0)); // 0.954
  HELIB_STATS_UPDATE("sigma_2_5", double(std_devs <= 2.5)); // 0.988
  HELIB_STATS_UPDATE("sigma_3_0", double(std_devs <= 3.0)); // 0.997, 1 in 370
  HELIB_STATS_UPDATE("sigma_3_5",
                     double(std_devs <= 3.5)); // 0.999535, 1 in 2149
  HELIB_STATS_UPDATE("sigma_4_0",
                     double(std_devs <= 4.0)); // 0.999937, 1 in 15787

  // compute sample variance, and scale by the variance we expect
  HELIB_STATS_UPDATE("sigma_calc",
                     fsquare(ran_pwrfl) / fsquare(p2ePrime * sigma));

  // save the scaled value for application of other tests
  HELIB_STATS_SAVE("v_values", ran_pwrfl / (p2ePrime * sigma));
}

#endif

#if 0
void fhe_stats_print(long iter, const Context& context)
{
   long phim = context.getPhiM();

   std::cerr << "||||| recryption stats ||||\n";
   std::cerr << "**** averages ****\n";
   std::cerr << "=== critical_value=" << (fhe_stats_cv_sum/iter) << "\n";
   std::cerr << "=== |x|/bound=" << (fhe_stats_x_sum/iter) << "\n";
   std::cerr << "=== |x%q|/bound=" << (fhe_stats_xmod_sum/iter) << "\n";
   std::cerr << "=== |u|/bound=" << (fhe_stats_u_sum/iter) << "\n";
   std::cerr << "=== |v|/bound=" << (fhe_stats_v_sum/iter) << "\n";
   std::cerr << "**** maxima ****\n";
   std::cerr << "=== critical_value=" << (fhe_stats_cv_max) << "\n";
   std::cerr << "=== |x|/bound=" << (fhe_stats_x_max) << "\n";
   std::cerr << "=== |x%q|/bound=" << (fhe_stats_xmod_max) << "\n";
   std::cerr << "=== |u|/bound=" << (fhe_stats_u_max) << "\n";
   std::cerr << "=== |v|/bound=" << (fhe_stats_v_max) << "\n";
   std::cerr << "**** theoretical bounds ***\n";
   std::cerr << "=== single-max=" << (sqrt(2.0*log(phim))/context.scale) << "\n";
   std::cerr << "=== global-max=" << (sqrt(2.0*(log(iter)+log(phim)))/context.scale) << "\n";


}
#endif

} // namespace helib
