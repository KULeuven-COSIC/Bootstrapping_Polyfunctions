/* Copyright (C) 2012-2021 IBM Corp.
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
// Notice: this file was modified from HElib
#ifndef HELIB_CTXT_H
#define HELIB_CTXT_H
/**
 * @file Ctxt.h
 * @brief Declarations of a BGV-type ciphertext and key-switching matrices
 *
 * A ciphertext is a std::vector of "ciphertext parts", each part consists of
 * a polynomial (element of polynomial ring R_Q) and a "handle" describing
 * the secret-key polynomial that this part multiplies during decryption.
 * For example:
 * + A "canonical" ciphertext has two parts, the first part multiplies 1
 *   and the second multiplies the "base" secret key s.
 * + When you multiply two canonical ciphertexts you get a 3-part
 *   ciphertext, with parts corresponding to 1, s, and s^2.
 * + When you apply automorphism X->X^t to a generic ciphertext, then
 *   - the part corresponding to 1 still remains wrt 1
 *   - every other part corresponding to some s' will now be corresponding
 *     to the polynomial s'(X^t) mod Phi_m(X)
 *
 * This type of representation lets you in principle add ciphertexts that
 * are defined with respect to different keys:
 * + For parts of the two ciphertexts that point to the same secret-key
 *   polynomial, you just add the two Double-CRT polynomials
 * + Parts in one ciphertext that do not have counter-part in the other
 *   ciphertext will just be included in the result intact.
 * For example, you have the ciphertexts
 *    C1 = (a relative to 1, b relative to s)
 *    C2 = (u relative to 1, v relative to s(X^3))
 * Then their sum will be
 *    C1+C2 = (a+u relative to 1, b relative to s, v relative to s(X^3))
 *
 * Similarly, in principle you can also multiply arbitrary ciphertexts, even
 * ones that are defined with respect to different keys, and the result will
 * be defined with respect to the tensor product of the two keys.
 *
 * The current implementation is more restrictive, however. It requires that
 * a ciphertext has one part wrt 1, that for every r>=1 there is at most one
 * part wrt to s^r(X^t) (for some t), and that the r's are consecutive. For
 * example you cannot have parts wrt (1,s,s^3) without having a part wrt s^2.
 *
 * It follows that you can only add/multiply ciphertexts if one of the two
 * lists of handles is a prefix of the other. For example, one can add a
 * ciphertext wrt (1,s(X^2)) to another wrt (1,s(X^2),s^2(X^2)), but not
 * to another ciphertext wrt (1,s).
 **/

#include <helib/DoubleCRT.h>
#include <helib/EncodedPtxt.h>
#include <helib/apiAttributes.h>

#include <cfloat> // DBL_MAX

namespace helib {
struct CKKS;
struct BGV;

template <typename Scheme>
class Ptxt;

class KeySwitch;
class PubKey;
class SecKey;

class PtxtArray;

/**
 * @class SKHandle
 * @brief A handle, describing the secret-key element that "matches" a part, of
 * the form s^r(X^t).
 **/
class SKHandle
{
  long powerOfS, powerOfX, secretKeyID;

public:
  friend class Ctxt;

  SKHandle(long newPowerOfS = 0, long newPowerOfX = 1, long newSecretKeyID = 0)
  {
    powerOfS = newPowerOfS;
    powerOfX = newPowerOfX;
    secretKeyID = newSecretKeyID;
  }

  //! @brief Set powerOfS=powerOfX=1
  void setBase(long newSecretKeyID = -1)
  {
    powerOfS = 1;
    powerOfX = 1;
    if (newSecretKeyID >= 0)
      secretKeyID = newSecretKeyID;
  }

  //! @brief Is powerOfS==powerOfX==1?
  bool isBase(long ofKeyID = 0) const
  {
    // If ofKeyID<0, only check that this is base of some key,
    // otherwise check that this is base of the given key
    return powerOfS == 1 && powerOfX == 1 &&
           (ofKeyID < 0 || secretKeyID == ofKeyID);
  }

  //! @brief Set powerOfS=0, powerOfX=1
  void setOne(long newSecretKeyID = -1)
  {
    powerOfS = 0;
    powerOfX = 1;
    if (newSecretKeyID >= 0)
      secretKeyID = newSecretKeyID;
  }

  //! @brief Is powerOfS==0?
  bool isOne() const { return powerOfS == 0; }

  bool operator==(const SKHandle& other) const
  {
    if (powerOfS == 0 && other.powerOfS == 0)
      return true;
    else
      return powerOfS == other.powerOfS && powerOfX == other.powerOfX &&
             secretKeyID == other.secretKeyID;
  }

  bool operator!=(const SKHandle& other) const { return !(*this == other); }

  /* access functions */

  long getPowerOfS() const { return powerOfS; }
  long getPowerOfX() const { return powerOfX; }
  long getSecretKeyID() const { return secretKeyID; }

  /**
   * @brief Computes the "product" of two handles.
   *
   * The key-ID's and powers of X must match, else an error state arises,
   * which is represented using a key-ID of -1 and returning false. Also,
   * note that inputs may alias outputs.
   *
   * To determine if the resulting handle can be re-linearized using
   * some key-switching matrices from the public key, use the method
   * pubKey.haveKeySWmatrix(handle,handle.secretKeyID), from the class
   * PubKey in keys.h
   */
  bool mul(const SKHandle& a, const SKHandle& b)
  {
    // If either of the inputs is one, the output equals to the other input
    if (a.isOne()) {
      *this = b;
      return (b.secretKeyID >= 0);
    }
    if (b.isOne()) {
      *this = a;
      return (a.secretKeyID >= 0);
    }

    if (a.secretKeyID == -1 || b.secretKeyID == -1) {
      secretKeyID = -1; // -1 will be used to indicate an "error state"
      return false;
    }

    if (a.secretKeyID != b.secretKeyID) {
      secretKeyID = -1;
      return false;
    }

    if (a.powerOfX != b.powerOfX) {
      secretKeyID = -1;
      return false;
    }

    secretKeyID = a.secretKeyID;
    powerOfX = a.powerOfX;
    powerOfS = a.powerOfS + b.powerOfS;
    return true;
  }

  friend std::istream& operator>>(std::istream& s, SKHandle& handle);

  // Raw IO

  /**
   * @brief Write out the `SKHandle` object in binary format.
   * @param str Output `std::ostream`.
   **/
  void writeTo(std::ostream& str) const;

  /**
   * @brief Read from the stream the serialized `SKHandle` object in binary
   * format.
   * @param str Input `std::istream`.
   * @return The deserialized `SKHandle` object.
   **/
  static SKHandle readFrom(std::istream& str);

  /**
   * @brief Write out the secret key handle (`SKHandle`) object to the output
   * stream using JSON format.
   * @param str Output `std::ostream`.
   **/
  void writeToJSON(std::ostream& str) const;

  /**
   * @brief Write out the secret key handle (`SKHandle`) object to a
   * `JsonWrapper`.
   * @return The `JsonWrapper`.
   **/
  JsonWrapper writeToJSON() const;

  /**
   * @brief Read from the stream the serialized secret key handle (`SKHandle`)
   * object using JSON format.
   * @param str Input `std::istream`.
   * @return The deserialized `SKHandle` object.
   **/
  static SKHandle readFromJSON(std::istream& str);

  /**
   * @brief Read from the `JsonWrapper` the serialized secret key handle
   * (`SKHandle`) object.
   * @param j The `JsonWrapper` containing the serialized `SKHandle` object.
   * @return The deserialized `SKHandle` object.
   **/
  static SKHandle readFromJSON(const JsonWrapper& j);

  /**
   * @brief Read from the stream the serialized secret key handle (`SKHandle`)
   * object using JSON format.
   * @param str Input `std::istream`.
   **/
  void readJSON(std::istream& str);

  /**
   * @brief Read from the `JsonWrapper` the serialized secret key handle
   * (`SKHandle`) object.
   * @param j The `JsonWrapper` containing the serialized `SKHandle` object.
   **/
  void readJSON(const JsonWrapper& j);
};

inline std::ostream& operator<<(std::ostream& s, const SKHandle& handle)
{
  handle.writeToJSON(s);
  return s;
}

/**
 * @class CtxtPart
 * @brief One entry in a ciphertext std::vector
 *
 * A ciphertext part consists of a polynomial (element of the ring R_Q)
 * and a handle to the corresponding secret-key polynomial.
 **/
class CtxtPart : public DoubleCRT
{
public:
  //! @brief The handle is a public data member
  SKHandle skHandle; // The secret-key polynomial corresponding to this part

  bool operator==(const CtxtPart& other) const;
  bool operator!=(const CtxtPart& other) const { return !(*this == other); }

  // Constructors

  CtxtPart(const Context& _context, const IndexSet& s) : DoubleCRT(_context, s)
  {
    skHandle.setOne();
  }

  CtxtPart(const Context& _context,
           const IndexSet& s,
           const SKHandle& otherHandle) :
      DoubleCRT(_context, s), skHandle(otherHandle)
  {}

  // Copy constructors from the base class
  explicit CtxtPart(const DoubleCRT& other) : DoubleCRT(other)
  {
    skHandle.setOne();
  }

  CtxtPart(const DoubleCRT& other, const SKHandle& otherHandle) :
      DoubleCRT(other), skHandle(otherHandle)
  {}

  /**
   * @brief Write out the `CtxtPart` object in binary format.
   * @param str Output `std::ostream`.
   **/
  void writeTo(std::ostream& str) const;

  /**
   * @brief Read from the stream the serialized `CtxtPart` object in binary
   * format.
   * @param str Input `std::istream`.
   * @return The deserialized `CtxtPart` object.
   **/
  static CtxtPart readFrom(std::istream& str, const Context& context);

  /**
   * @brief In-place read from the stream the serialized `CtxtPart` object in
   * binary format.
   * @param str Input `std::istream`.
   **/
  void read(std::istream& str);

  /**
   * @brief Write out the ciphertext part (`CtxtPart`) object to the output
   * stream using JSON format.
   * @param str Output `std::ostream`.
   **/
  void writeToJSON(std::ostream& str) const;

  /**
   * @brief Write out the ciphertext part (`CtxtPart`) object to a
   * `JsonWrapper`.
   * @return The `JsonWrapper`.
   **/
  JsonWrapper writeToJSON() const;

  /**
   * @brief Read from the stream the serialized ciphertext part (`CtxtPart`)
   * object using JSON format.
   * @param str Input `std::istream`.
   * @param context The `Context` to be used.
   * @return The deserialized `CtxtPart` object.
   **/
  static CtxtPart readFromJSON(std::istream& str, const Context& context);

  /**
   * @brief Read from the `JsonWrapper` the serialized ciphertext part
   *(`CtxtPart`) object.
   * @param j The `JsonWrapper` containing the serialized `CtxtPart` object.
   * @param context The `Context` to be used.
   * @return The deserialized `CtxtPart` object.
   **/
  static CtxtPart readFromJSON(const JsonWrapper& j, const Context& context);

  /**
   * @brief Read from the stream the serialized ciphertext part (`CtxtPart`)
   * object using JSON format.
   * @param str Input `std::istream`.
   **/
  void readJSON(std::istream& str);

  /**
   * @brief Read from the `JsonWrapper` the serialized ciphertext part
   *(`CtxtPart`) object.
   * @param j The `JsonWrapper` containing the serialized `SKHandle` object.
   **/
  void readJSON(const JsonWrapper& j);
};

std::istream& operator>>(std::istream& s, CtxtPart& p);
std::ostream& operator<<(std::ostream& s, const CtxtPart& p);

//! \cond FALSE (make doxygen ignore this code)
struct ZeroCtxtLike_type
{}; // used to select a constructor

const ZeroCtxtLike_type ZeroCtxtLike = ZeroCtxtLike_type();
//! \endcond

/**
 * @class Ctxt
 * @brief A Ctxt object holds a single ciphertext
 *
 * The class Ctxt includes a std::vector<CtxtPart>: For a Ctxt c, c[i] is the
 * i'th ciphertext part, which can be used also as a DoubleCRT object (since
 * CtxtPart is derived from DoubleCRT). By convention, c[0], the first CtxtPart
 * object in the std::vector, has skHndl that points to 1 (i.e., it is just
 * added in upon decryption, without being multiplied by anything).  We
 * maintain the invariance that all the parts of a ciphertext are defined
 * relative to the same set of primes.
 *
 * A ciphertext contains also pointers to the general parameters of this FHE
 * instance and the public key, and a high-probability bound on the noise
 * magnitude (kept in the noiseBound data member). The noise bound is a bound
 * on the l-infinity norm of the canonical embedding of the noise polynomial,
 * namely its evaluation in roots of the ring polynomial (which are the complex
 * primitive roots of unity).  The noise bound is added on addition, multiplied
 * on multiplications, remains unchanged for automorphism, and is roughly
 * scaled down by mod-switching with some added factor, and similarly scaled up
 * by key-switching with some added factor.
 *
 **/
class Ctxt
{
  friend class PubKey;
  friend class SecKey;
  friend class BasicAutomorphPrecon;

  const Context& context;      // points to the parameters of this FHE instance
  const PubKey& pubKey;        // points to the public encryption key;
  std::vector<CtxtPart> parts; // the ciphertext parts
  IndexSet primeSet; // the primes relative to which the parts are defined
  long ptxtSpace;    // plaintext space for this ciphertext (either p or p^r)

  // a high-probability bound on the noise magnitude
  NTL::xdouble noiseBound;

  long intFactor; // an integer factor to divide by on decryption (for BGV)
  NTL::xdouble ratFactor; // rational factor to divide on decryption (for CKKS)
  NTL::xdouble ptxtMag;   // bound on the plaintext size (for CKKS)

  // Create a tensor product of c1,c2. It is assumed that *this,c1,c2
  // are defined relative to the same set of primes and plaintext space,
  // and that *this DOES NOT point to the same object as c1,c2
  void tensorProduct(const Ctxt& c1, const Ctxt& c2);

  // Add/subtract a ciphertext part to/from a ciphertext. These are private
  // methods, they cannot update the noiseBound so they must be called
  // from a procedure that will eventually update that estimate.
  Ctxt& operator-=(const CtxtPart& part)
  {
    subPart(part);
    return *this;
  }

  Ctxt& operator+=(const CtxtPart& part)
  {
    addPart(part);
    return *this;
  }

  // NOTE: the matchPrimeSets business in the following routines
  // is DEPRECATED.  The requirement is that the prime set of part
  // must contain the prime set of *this.
  // If not, an exception is raised.
  // Also, if matchPrimeSets == true and the prime set of *this does
  // not contain the prime set of part, an exception is also raised.
  // Also, note that the prime set of *this will not change.

  // Procedural versions with additional parameter
  void subPart(const CtxtPart& part, bool matchPrimeSet = false)
  {
    subPart(part, part.skHandle, matchPrimeSet);
  }

  void addPart(const CtxtPart& part, bool matchPrimeSet = false)
  {
    addPart(part, part.skHandle, matchPrimeSet);
  }

  void subPart(const DoubleCRT& part,
               const SKHandle& handle,
               bool matchPrimeSet = false)
  {
    addPart(part, handle, matchPrimeSet, true);
  }

  void addPart(const DoubleCRT& part,
               const SKHandle& handle,
               bool matchPrimeSet = false,
               bool negative = false);

  // convenient to avoid dealing with the deprecated matchPrimeSet
  // parameter
  void addSignedPart(const DoubleCRT& part,
                     const SKHandle& handle,
                     bool negative = false)
  {
    addPart(part, handle, false, negative);
  }

  // Takes as arguments a ciphertext-part p relative to s' and a key-switching
  // matrix W = W[s'->s], use W to switch p relative to (1,s), and add the
  // result to *this.
  void keySwitchPart(const CtxtPart& p, const KeySwitch& W);

  // internal procedure used in key-switching
  void keySwitchDigits(const KeySwitch& W, std::vector<DoubleCRT>& digits);

  long getPartIndexByHandle(const SKHandle& handle) const
  {
    for (size_t i = 0; i < parts.size(); i++)
      if (parts[i].skHandle == handle)
        return i;
    return -1;
  }

  // Sanity-check: Check that prime-set is "valid", i.e. that it
  // contains either all the special primes or none of them
  bool verifyPrimeSet() const;

  // A private assignment method that does not check equality of context or
  // public key, this is needed when we copy the pubEncrKey member between
  // different public keys.
  Ctxt& privateAssign(const Ctxt& other);

  // explicitly multiply intFactor by e, which should be
  // in the interval [0, ptxtSpace)
  void mulIntFactor(long e);

public:
  /**
   * @brief Class label to be added to JSON serialization as object type
   * information.
   */
  static constexpr std::string_view typeName = "Ctxt";

  // Default copy-constructor
  Ctxt(const Ctxt& other) = default;

  // VJS-FIXME: this was really a messy design choice to not
  // have ciphertext constructors that specify prime sets.
  // The default value of ctxtPrimes is kind of pointless.

  //__attribute__((deprecated))
  explicit Ctxt(const PubKey& newPubKey, long newPtxtSpace = 0); // constructor

  //__attribute__((deprecated))
  Ctxt(ZeroCtxtLike_type, const Ctxt& ctxt);
  // constructs a zero ciphertext with same public key and
  // plaintext space as ctxt

  //! Dummy encryption, just encodes the plaintext in a Ctxt object
  //! If provided, size should be a high-probability bound
  //! on the L-infty norm of the canonical embedding
  void DummyEncrypt(const NTL::ZZX& ptxt, double size = -1.0);

  Ctxt& operator=(const Ctxt& other)
  { // public assignment operator
    assertEq(&context,
             &other.context,
             "Cannot assign Ctxts with different context");
    assertEq(&pubKey,
             &other.pubKey,
             "Cannot assign Ctxts with different pubKey");
    return privateAssign(other);
  }

  bool operator==(const Ctxt& other) const { return equalsTo(other); }
  bool operator!=(const Ctxt& other) const { return !equalsTo(other); }

  // a procedural variant with an additional parameter
  bool equalsTo(const Ctxt& other, bool comparePkeys = true) const;

  // Encryption and decryption are done by the friends [Pub|Sec]Key

  //! @name Ciphertext arithmetic
  ///@{
  void negate();

  // Add/subtract another ciphertext
  Ctxt& operator+=(const Ctxt& other)
  {
    addCtxt(other);
    return *this;
  }

  Ctxt& operator-=(const Ctxt& other)
  {
    addCtxt(other, true);
    return *this;
  }

  void addCtxt(const Ctxt& other, bool negative = false);

  // Multiply by another ciphertext
  void multLowLvl(const Ctxt& other, bool destructive = false);

  // This is a high-level mul with relinearization
  Ctxt& operator*=(const Ctxt& other)
  {
    multiplyBy(other);
    return *this;
  }

  void automorph(long k); // Apply automorphism F(X) -> F(X^k) (gcd(k,m)=1)
  Ctxt& operator>>=(long k)
  {
    automorph(k);
    return *this;
  }

  void complexConj(); // Complex conjugate, same as automorph(m-1)

  //! @brief automorphism with re-linearization
  void smartAutomorph(long k);
  // Apply F(X)->F(X^k) followed by re-linearization. The automorphism is
  // possibly evaluated via a sequence of steps, to ensure that we can
  // re-linearize the result of every step.

  //! @brief applies the automorphism p^j using smartAutomorphism
  void frobeniusAutomorph(long j);

  /**
   * @brief Times equals operator with a `ZZX`.
   * @param poly Element by which to multiply.
   * @return Reference to `*this` post multiplication.
   * @deprecated This function is deprecated in favor of a new
   * `EncodedPtxt`-based API.\n
   * Please use `Ctxt::operator*=(const EncodedPtxt& ptxt)` instead.
   **/
  [[deprecated(
      "Please use Ctxt::operator*=(const EncodedPtxt& ptxt) instead.")]] Ctxt&
  operator*=(const NTL::ZZX& poly);

  //! Add a constant polynomial.
  //! If provided, size should be a high-probability bound
  //! on the L-infty norm of the canonical embedding
  //! Otherwise, for the DoubleCRT variant, a bound based on the assumption
  //! that the coefficients are uniformly and independently distributed over
  //! [-ptxtSpace/2, ptxtSpace/2].
  //! For the other variants, explicit bounds are computed (if not CKKS).
  void addConstant(const DoubleCRT& dcrt, double size = -1.0);
  void addConstant(const NTL::ZZX& poly, double size = -1.0);

  /**
   * @brief Add a plaintext to this `Ctxt`.
   * @param ptxt Plaintext `Ptxt` object with which to add.
   **/
  template <typename Scheme>
  void addConstant(const Ptxt<Scheme>& ptxt, bool neg = false)
  {
    EncodedPtxt eptxt;
    ptxt.encode(eptxt);
    addConstant(eptxt, neg);
  }

  template <typename Scheme>
  Ctxt& operator+=(const Ptxt<Scheme>& ptxt)
  {
    addConstant(ptxt);
    return *this;
  }

  template <typename Scheme>
  Ctxt& operator-=(const Ptxt<Scheme>& ptxt)
  {
    addConstant(ptxt, true);
    return *this;
  }

  //! add a rational number in the form a/b, a,b are long
  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator+=(double ptxt)` instead.
   **/
  [[deprecated("Please use Ctxt::operator+=(double ptxt) instead.")]] void
      addConstantCKKS(std::pair</*numerator=*/long, /*denominator=*/long>);
  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator+=(double ptxt)` instead.
   **/
  [[deprecated("Please use Ctxt::operator+=(double ptxt) instead.")]] void
  addConstantCKKS(double x)
  { // FIXME: not enough precision when x is large
    // This function is deprecated.
    //    addConstantCKKS(
    //        rationalApprox(x, /*denomBound=*/1L <<
    //        getContext().getAlMod().getR()));
    auto p =
        rationalApprox(x, /*denomBound=*/1L << getContext().getAlMod().getR());
    *this += double(p.first) / p.second;
  }

  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator+=(const FatEncodedPtxt& ptxt)` instead.
   **/
  // [[deprecated]]
  void addConstantCKKS(const DoubleCRT& dcrt,
                       NTL::xdouble size = NTL::xdouble(-1.0),
                       NTL::xdouble factor = NTL::xdouble(-1.0));

  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator+=(const EncodedPtxt& ptxt)` instead.
   **/
  [[deprecated(
      "Please use Ctxt::operator+=(const EncodedPtxt& ptxt) instead.")]] void
  addConstantCKKS(const NTL::ZZX& poly,
                  NTL::xdouble size = NTL::xdouble(-1.0),
                  NTL::xdouble factor = NTL::xdouble(-1.0));

  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator+=(const PtxtArray& ptxt)` instead.
   **/
  [[deprecated(
      "Please use Ctxt::operator+=(const PtxtArray& ptxt) instead.")]] void
  addConstantCKKS(const std::vector<std::complex<double>>& ptxt);

  /**
   * @brief Add a `CKKS` plaintext to this `Ctxt`.
   * @param ptxt Plaintext `Ptxt` object with which to add.
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator+=(const PtxtArray& ptxt)` instead.
   **/
  [[deprecated(
      "Please use Ctxt::operator+=(const Ptxt<Scheme>& ptxt) instead.")]] void
  addConstantCKKS(const Ptxt<CKKS>& ptxt);
  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator+=(const NTL::ZZ& ptxt)` instead.
   **/
  [[deprecated(
      "Please use Ctxt::operator+=(const NTL::ZZ& ptxt) instead.")]] void
  addConstantCKKS(const NTL::ZZ& c);

  //! Multiply-by-constant.
  //! If the size is not given, for the DCRT variant, we use
  //! a high probability bound assuming "random" coefficients
  //! mod ptxtSpace, while for the other variants, we use
  //! explicitly computed bounds (if not CKKS).
  void multByConstant(const DoubleCRT& dcrt, double size = -1.0);

  void multByConstant(const NTL::ZZX& poly, double size = -1.0);
  void multByConstant(const zzX& poly, double size = -1.0);

  //=========== new multByConstant interface =========

  /**
   * @brief Multiply a `Ctxt` with a specified plaintext constant.
   * @param ptxt The constant to multiply as a `PtxtArray` object.
   **/
  void multByConstant(const PtxtArray& ptxt);
  /**
   * @brief Multiply a `Ctxt` with a specified plaintext constant.
   * @param ptxt The constant to multiply as an `EncodedPtxt` object.
   * @note `EncodedPtxt` is a plaintext object containing `NTL::ZZX` data.
   **/
  void multByConstant(const EncodedPtxt& ptxt);
  /**
   * @brief Multiply a `Ctxt` with a specified plaintext constant.
   * @param ptxt The constant to multiply as a `FatEncodedPtxt` object.
   * @note `FatEncodedPtxt` is a plaintext object containing
   * `helib::DoubleCRT` data.
   **/
  void multByConstant(const FatEncodedPtxt& ptxt);

  /**
   * @brief Multiply a `Ctxt` with an `NTL::ZZ` scalar.
   * @param ptxt Scalar to multiply.
   **/
  void multByConstant(const NTL::ZZ& ptxt);
  /**
   * @brief Multiply a `Ctxt` with a `long` scalar.
   * @param ptxt Scalar to multiply.
   **/
  void multByConstant(long ptxt);
  /**
   * @brief Multiply a `Ctxt` with a `double` scalar.
   * @param ptxt Scalar to multiply.
   **/
  void multByConstant(double ptxt);
  /**
   * @brief Multiply a `Ctxt` with an `NTL::xdouble` scalar.
   * @param ptxt Scalar to multiply.
   **/
  void multByConstant(NTL::xdouble ptxt);

  /**
   * @brief Times equals operator with a plaintext constant.
   * @param ptxt Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   **/
  Ctxt& operator*=(const PtxtArray& ptxt)
  {
    multByConstant(ptxt);
    return *this;
  }

  /**
   * @brief Times equals operator with a plaintext constant.
   * @param ptxt Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   * @note `EncodedPtxt` is a plaintext object containing `NTL::ZZX` data.
   **/
  Ctxt& operator*=(const EncodedPtxt& ptxt)
  {
    multByConstant(ptxt);
    return *this;
  }

  /**
   * @brief Times equals operator with a plaintext constant.
   * @param ptxt Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   * @note `FatEncodedPtxt` is a plaintext object containing
   * `helib::DoubleCRT` data.
   **/
  Ctxt& operator*=(const FatEncodedPtxt& ptxt)
  {
    multByConstant(ptxt);
    return *this;
  }

  /**
   * @brief Times equals operator with an `NTL::ZZ` scalar.
   * @param ptxt Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   **/
  Ctxt& operator*=(const NTL::ZZ& ptxt)
  {
    multByConstant(ptxt);
    return *this;
  }

  /**
   * @brief Times equals operator with a `long` scalar.
   * @param ptxt Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   **/
  Ctxt& operator*=(long ptxt)
  {
    multByConstant(ptxt);
    return *this;
  }

  /**
   * @brief Times equals operator with a `double` scalar.
   * @param ptxt Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   **/
  Ctxt& operator*=(double ptxt)
  {
    multByConstant(ptxt);
    return *this;
  }

  /**
   * @brief Times equals operator with an `NTL::xdouble` scalar.
   * @param ptxt Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   **/
  Ctxt& operator*=(NTL::xdouble ptxt)
  {
    multByConstant(ptxt);
    return *this;
  }

private: // impl only
  void multByConstant(const FatEncodedPtxt_BGV& ptxt);
  void multByConstant(const FatEncodedPtxt_CKKS& ptxt);

public:
  //=========== new addConstant interface ============

  /**
   * @brief Add to a `Ctxt` a specified plaintext constant.
   * @param ptxt The constant to add as a `PtxtArray` object.
   * @param neg Flag to specify if the constant is negative. Default is
   * `false`.
   **/
  void addConstant(const PtxtArray& ptxt, bool neg = false);
  /**
   * @brief Add to a `Ctxt` a specified plaintext constant.
   * @param ptxt The constant to add as an `EncodedPtxt` object.
   * @param neg Flag to specify if the constant is negative. Default is
   * `false`.
   * @note `EncodedPtxt` is a plaintext object containing `NTL::ZZX` data.
   **/
  void addConstant(const EncodedPtxt& ptxt, bool neg = false);
  /**
   * @brief Add to a `Ctxt` a specified plaintext constant.
   * @param ptxt The constant to add as a `FatEncodedPtxt` object.
   * @param neg Flag to specify if the constant is negative. Default is
   * `false`.
   * @note `FatEncodedPtxt` is a plaintext object containing
   * `helib::DoubleCRT` data.
   **/
  void addConstant(const FatEncodedPtxt& ptxt, bool neg = false);

  /**
   * @brief Add to a `Ctxt` an `NTL::ZZ` scalar.
   * @param ptxt The scalar to add.
   * @param neg Flag to specify if the constant is negative. Default is
   * `false`.
   **/
  void addConstant(const NTL::ZZ& ptxt, bool neg = false);
  /**
   * @brief Add to a `Ctxt` a `long` scalar.
   * @param ptxt The scalar to add.
   * @param neg Flag to specify if the constant is negative. Default is
   * `false`.
   **/
  void addConstant(long ptxt, bool neg = false);
  /**
   * @brief Add to a `Ctxt` a `double` scalar.
   * @param ptxt The scalar to add.
   * @param neg Flag to specify if the constant is negative. Default is
   * `false`.
   **/
  void addConstant(double ptxt, bool neg = false);
  /**
   * @brief Add to a `Ctxt` an `NTL::xdouble` scalar.
   * @param ptxt The scalar to add.
   * @param neg Flag to specify if the constant is negative. Default is
   * `false`.
   **/
  void addConstant(NTL::xdouble ptxt, bool neg = false);

  /**
   * @brief Plus equals operator with plaintext constant.
   * @param ptxt Right hand side of addition.
   * @return Reference to `*this` post addition.
   **/
  Ctxt& operator+=(const PtxtArray& ptxt)
  {
    addConstant(ptxt);
    return *this;
  }

  /**
   * @brief Plus equals operator with plaintext constant.
   * @param ptxt Right hand side of addition.
   * @return Reference to `*this` post addition.
   * @note `EncodedPtxt` is a plaintext object containing `NTL::ZZX` data.
   **/
  Ctxt& operator+=(const EncodedPtxt& ptxt)
  {
    addConstant(ptxt);
    return *this;
  }

  /**
   * @brief Plus equals operator with plaintext constant.
   * @param ptxt Right hand side of addition.
   * @return Reference to `*this` post addition.
   * @note `FatEncodedPtxt` is a plaintext object containing
   * `helib::DoubleCRT` data.
   **/
  Ctxt& operator+=(const FatEncodedPtxt& ptxt)
  {
    addConstant(ptxt);
    return *this;
  }

  /**
   * @brief Plus equals operator with an `NTL::ZZ` scalar.
   * @param ptxt Right hand side of addition.
   * @return Reference to `*this` post addition.
   **/
  Ctxt& operator+=(const NTL::ZZ& ptxt)
  {
    addConstant(ptxt);
    return *this;
  }

  /**
   * @brief Plus equals operator with a `long` scalar.
   * @param ptxt Right hand side of addition.
   * @return Reference to `*this` post addition.
   **/
  Ctxt& operator+=(long ptxt)
  {
    addConstant(ptxt);
    return *this;
  }

  /**
   * @brief Plus equals operator with a `double` scalar.
   * @param ptxt Right hand side of addition.
   * @return Reference to `*this` post addition.
   **/
  Ctxt& operator+=(double ptxt)
  {
    addConstant(ptxt);
    return *this;
  }

  /**
   * @brief Plus equals operator with an `NTL::xdouble` scalar.
   * @param ptxt Right hand side of addition.
   * @return Reference to `*this` post addition.
   **/
  Ctxt& operator+=(NTL::xdouble ptxt)
  {
    addConstant(ptxt);
    return *this;
  }

  /**
   * @brief Minus equals operator with plaintext constant.
   * @param ptxt Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   **/
  Ctxt& operator-=(const PtxtArray& ptxt)
  {
    addConstant(ptxt, true);
    return *this;
  }

  /**
   * @brief Minus equals operator with plaintext constant.
   * @param ptxt Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   * @note `EncodedPtxt` is a plaintext object containing `NTL::ZZX` data.
   **/
  Ctxt& operator-=(const EncodedPtxt& ptxt)
  {
    addConstant(ptxt, true);
    return *this;
  }

  /**
   * @brief Minus equals operator with plaintext constant.
   * @param ptxt Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   * @note `FatEncodedPtxt` is a plaintext object containing
   * `helib::DoubleCRT` data.
   **/
  Ctxt& operator-=(const FatEncodedPtxt& ptxt)
  {
    addConstant(ptxt, true);
    return *this;
  }

  /**
   * @brief Minus equals operator with an `NTL::ZZ` scalar.
   * @param ptxt Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   **/
  Ctxt& operator-=(const NTL::ZZ& ptxt)
  {
    addConstant(ptxt, true);
    return *this;
  }

  /**
   * @brief Minus equals operator with a `long` scalar.
   * @param ptxt Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   **/
  Ctxt& operator-=(long ptxt)
  {
    addConstant(ptxt, true);
    return *this;
  }

  /**
   * @brief Minus equals operator with a `double` scalar.
   * @param ptxt Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   **/
  Ctxt& operator-=(double ptxt)
  {
    addConstant(ptxt, true);
    return *this;
  }

  /**
   * @brief Minus equals operator with an `NTL::xdouble` scalar.
   * @param ptxt Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   **/
  Ctxt& operator-=(NTL::xdouble ptxt)
  {
    addConstant(ptxt, true);
    return *this;
  }

private: // impl only
  void addConstant(const FatEncodedPtxt_BGV& ptxt, bool neg = false);
  void addConstant(const EncodedPtxt_BGV& ptxt, bool neg = false);

  void addConstant(const FatEncodedPtxt_CKKS& ptxt, bool neg = false);

public:
  //==================================================

  /**
   * @brief Multiply a plaintext to this `Ctxt`.
   * @param ptxt Plaintext `Ptxt` object with which to multiply.
   **/
  template <typename Scheme>
  void multByConstant(const Ptxt<Scheme>& ptxt)
  {
    EncodedPtxt eptxt;
    ptxt.encode(eptxt);
    multByConstant(eptxt);
  }

  template <typename Scheme>
  Ctxt& operator*=(const Ptxt<Scheme>& ptxt)
  {
    multByConstant(ptxt);
    return *this;
  }

  //! multiply by a rational number or floating point
  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator*=(double ptxt)` instead.
   **/
  [[deprecated("Please use Ctxt::operator*=(double ptxt) instead.")]] void
  multByConstantCKKS(double x)
  {
    if (this->isEmpty())
      return;

    if (x == 1)
      return;

    if (x == 0) {
      clear();
      return;
    }

    double size = std::abs(x);
    ptxtMag *= size;
    ratFactor /= size;
    if (x < 0)
      this->negate();
  }

  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator*=(double ptxt)` instead.
   **/
  [[deprecated("Please use Ctxt::operator*=(double ptxt) instead.")]] void
  multByConstantCKKS(std::pair<long, long> num) // rational number
  {
    // This function is deprecated.
    // multByConstantCKKS(double(num.first) / num.second);
    *this *= double(num.first) / num.second;
  }

  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator*=(const FatEncodedPtxt& ptxt)` instead.
   **/
  // [[deprecated]]
  void multByConstantCKKS(const DoubleCRT& dcrt,
                          NTL::xdouble size = NTL::xdouble(-1.0),
                          NTL::xdouble factor = NTL::xdouble(-1.0),
                          double roundingErr = -1.0);

  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator*=(const EncodedPtxt& ptxt)` instead.
   **/
  [[deprecated(
      "Please use Ctxt::operator*=(const EncodedPtxt& ptxt) instead.")]] void
  multByConstantCKKS(const NTL::ZZX& poly,
                     NTL::xdouble size = NTL::xdouble(-1.0),
                     NTL::xdouble factor = NTL::xdouble(-1.0),
                     double roundingErr = -1.0)
  {
    DoubleCRT dcrt(poly, context, primeSet);
    multByConstantCKKS(dcrt, size, factor, roundingErr);
  }

  // TODO: Ptxt: refactor into single function for getting the
  // std::vector<cx_double> repr of slots
  /**
   * @brief Multiply a `CKKS` plaintext to this `Ctxt`.
   * @param ptxt Plaintext `Ptxt` object polynomial with which to multiply.
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator*=(const PtxtArray& ptxt)` instead.
   **/
  [[deprecated(
      "Please use Ctxt::operator*=(const Ptxt<Scheme>& ptxt) instead.")]] void
  multByConstantCKKS(const Ptxt<CKKS>& ptxt);
  /**
   * @deprecated This function is deprecated in favor of a newer API.
   * Please use `Ctxt::operator*=(const PtxtArray& ptxt)` instead.
   **/
  [[deprecated(
      "Please use Ctxt::operator*=(const PtxtArray& ptxt) instead.")]] void
  multByConstantCKKS(const std::vector<std::complex<double>>& ptxt);

  //! Convenience method: XOR and nXOR with arbitrary plaintext space:
  //! a xor b = a+b-2ab = a + (1-2a)*b,
  //! a nxor b = 1-a-b+2ab = (b-1)(2a-1)+a
  void xorConstant(const DoubleCRT& poly, UNUSED double size = -1.0)
  {
    DoubleCRT tmp = poly;
    tmp *= -2;
    tmp += 1; // tmp = 1-2*poly
    multByConstant(tmp);
    addConstant(poly);
  }

  void xorConstant(const NTL::ZZX& poly, double size = -1.0)
  {
    xorConstant(DoubleCRT(poly, context, primeSet), size);
  }

  void nxorConstant(const DoubleCRT& poly, UNUSED double size = -1.0)
  {
    DoubleCRT tmp = poly;
    tmp *= 2;
    tmp -= 1;                      // 2a-1
    addConstant(NTL::to_ZZX(-1L)); // b11
    multByConstant(tmp);           // (b-1)(2a-1)
    addConstant(poly);             // (b-1)(2a-1)+a = 1-a-b+2ab
  }

  void nxorConstant(const NTL::ZZX& poly, double size = -1.0)
  {
    nxorConstant(DoubleCRT(poly, context, primeSet), size);
  }

  //! Divide a ciphertext by p, for plaintext space p^r, r>1. It is assumed
  //! that the ciphertext encrypts a polynomial which is zero mod p. If this
  //! is not the case then the result will not be a valid ciphertext anymore.
  //! As a side-effect, the plaintext space is reduced from p^r to p^{r-1}.
  void divideByP();

  //! Multiply ciphertext by p^e, for plaintext space p^r. This also has
  //! the side-effect of increasing the plaintext space to p^{r+e}.
  void multByP(long e = 1)
  {
    long p2e = NTL::power_long(context.getP(), e);
    ptxtSpace *= p2e;
    multByConstant(NTL::to_ZZ(p2e));
  }

  // For backward compatibility
  void divideBy2();
  void extractBits(std::vector<Ctxt>& bits, long nBits2extract = 0);

  // Higher-level multiply routines
  void multiplyBy(const Ctxt& other);
  void customMultiplyBy(Ctxt& other, bool lazy);
  void multiplyBy2(const Ctxt& other1, const Ctxt& other2);
  void square() { multiplyBy(*this); }
  void cube() { multiplyBy2(*this, *this); }

  //! @brief raise ciphertext to some power
  void power(long e); // in polyEval.cpp

  ///@}

  //! @name Ciphertext maintenance
  ///@{

  //! Reduce plaintext space to a divisor of the original plaintext space
  void reducePtxtSpace(long newPtxtSpace);

  // This method can be used to increase the plaintext space, but the
  // high-order digits that you get this way are noise. Do not use it
  // unless you know what you are doing.
  void hackPtxtSpace(long newPtxtSpace) { ptxtSpace = newPtxtSpace; }

  void bumpNoiseBound(double factor) { noiseBound *= factor; }

  // CKKS adjustment to protect precision
  void relin_CKKS_adjust();

  void reLinearize(long keyIdx = 0);
  // key-switch to (1,s_i), s_i is the base key with index keyIdx

  Ctxt& cleanUp();
  // relinearize, then reduce, then drop special primes and small primes.

  // void reduce() const;

  //! @brief Add a high-noise encryption of the given constant
  void blindCtxt(const NTL::ZZX& poly);

  //! @brief Estimate the added noise
  NTL::xdouble modSwitchAddedNoiseBound() const;

  //! @brief Modulus-switching up (to a larger modulus).
  //! Must have primeSet <= s, and s must contain
  //! either all the special primes or none of them.
  void modUpToSet(const IndexSet& s);

  //! @brief Modulus-switching down (to a smaller modulus).
  //! mod-switch down to primeSet \intersect s, after this call we have
  //! primeSet<=s. s must contain either all special primes or none of them
  void modDownToSet(const IndexSet& s);

  //! @brief make the primeSet equal to newPrimeSet,
  //! via modUpToSet and modDownToSet
  void bringToSet(const IndexSet& s);

  // Finding the "natural" state of a ciphertext
  double naturalSize() const;       //! "natural size" is size before squaring
  IndexSet naturalPrimeSet() const; //! the corresponding primeSet

  //! @brief drop all smallPrimes and specialPrimes, adding ctxtPrimes
  //! as necessary to ensure that the scaled noise is above the
  //! modulus-switching added noise term.
  void dropSmallAndSpecialPrimes();

  //! @brief returns the *total* noise bound, which for CKKS
  //! is ptxtMag*ratFactor + noiseBound
  NTL::xdouble totalNoiseBound() const
  {
    if (isCKKS())
      return ptxtMag * ratFactor + noiseBound;
    else
      return noiseBound;
  }

  //! @brief for CKKS, returns a bound on the absolute error
  //! (which is noiseBound/ratFactor); for BGV, returns 0.
  double errorBound() const
  {
    if (isCKKS())
      return convert<double>(noiseBound / ratFactor);
    else
      return 0;
  }

  //! @brief returns the "capacity" of a ciphertext,
  //! which is the log2 of the ratio of the modulus to the
  //! *total* noise bound
  double capacity() const
  {
    return (logOfPrimeSet() -
            NTL::log(std::max(totalNoiseBound(), NTL::to_xdouble(1.0)))) /
           std::log(2.0);
  }

  //! @brief returns the log of the prime set
  double logOfPrimeSet() const { return context.logOfProduct(getPrimeSet()); }

  //! @brief the capacity in bits, returned as an integer
  long bitCapacity() const { return long(capacity()); }

  //! @brief Special-purpose modulus-switching for bootstrapping.
  //!
  //! Mod-switch to an externally-supplied modulus. The modulus need not be in
  //! the moduli-chain in the context, and does not even need to be a prime.
  //! The ciphertext *this is not affected, instead the result is returned in
  //! the zzParts std::vector, as a std::vector of ZZX'es.
  //! Returns an estimate for the scaled noise (not including the
  //! additive mod switching noise)
  double rawModSwitch(std::vector<NTL::ZZX>& zzParts, long toModulus) const;

  //! @brief compute the power X,X^2,...,X^n
  //  void computePowers(std::vector<Ctxt>& v, long nPowers) const;

  //! @brief Evaluate the cleartext poly on the encrypted ciphertext
  void evalPoly(const NTL::ZZX& poly);
  ///@}

  //! @name Utility methods
  ///@{

  void clear()
  { // set as an empty ciphertext
    parts.clear();
    primeSet = context.getCtxtPrimes();
    noiseBound = 0.0;
    intFactor = 1;
    ratFactor = ptxtMag = 1.0;
  }

  //! @brief Is this an empty ciphertext without any parts
  bool isEmpty() const { return (parts.size() == 0); }

  //! @brief A canonical ciphertext has (at most) handles pointing to (1,s)
  bool inCanonicalForm(long keyID = 0) const
  {
    if (parts.size() > 2)
      return false;
    if (parts.size() > 0 && !parts[0].skHandle.isOne())
      return false;
    if (parts.size() > 1 && !parts[1].skHandle.isBase(keyID))
      return false;
    return true;
  }

  //! @brief Would this ciphertext be decrypted without errors?
  bool isCorrect() const;

  const Context& getContext() const { return context; }
  const PubKey& getPubKey() const { return pubKey; }
  const IndexSet& getPrimeSet() const { return primeSet; }
  long getPtxtSpace() const { return ptxtSpace; }
  const NTL::xdouble& getNoiseBound() const { return noiseBound; }
  const NTL::xdouble& getRatFactor() const { return ratFactor; }
  const NTL::xdouble& getPtxtMag() const { return ptxtMag; }
  void setPtxtMag(const NTL::xdouble& z) { ptxtMag = z; }
  long getKeyID() const;

  bool isCKKS() const { return getContext().isCKKS(); }

  // Return r such that p^r = ptxtSpace
  long effectiveR() const;

  /**
   * @brief Returns log(noiseBound) - log(q)
   * @deprecated This is deprecated. Please use `Ctxt::capacity()` instead.
   **/
  [[deprecated("Please use Ctxt::capacity() instead.")]] double log_of_ratio()
      const
  {
    double logNoise =
        (getNoiseBound() <= 0.0) ? -DBL_MAX : log(getNoiseBound());
    double logMod =
        empty(getPrimeSet()) ? -DBL_MAX : context.logOfProduct(getPrimeSet());
    return logNoise - logMod;
  }
  ///@}
  friend std::istream& operator>>(std::istream& str, Ctxt& ctxt);
  friend std::ostream& operator<<(std::ostream& str, const Ctxt& ctxt);

  // Raw IO

  /**
   * @brief Write out the `Ctxt` object in binary format.
   * @param str Output `std::ostream`.
   **/
  void writeTo(std::ostream& str) const;

  /**
   * @brief Read from the stream the serialized `Ctxt` object in binary format.
   * @param str Input `std::istream`.
   * @return The deserialized `Ctxt` object.
   **/
  static Ctxt readFrom(std::istream& str, const PubKey& pubKey);

  /**
   * @brief In-place read from the stream the serialized `Ctxt` object in binary
   * format.
   * @param str Input `std::istream`.
   **/
  void read(std::istream& str);

  /**
   * @brief Write out the ciphertext (`Ctxt`) object to the output
   * stream using JSON format.
   * @param str Output `std::ostream`.
   **/
  void writeToJSON(std::ostream& str) const;

  /**
   * @brief Write out the ciphertext (`Ctxt`) object to a `JsonWrapper`.
   * @return The `JsonWrapper`.
   **/
  JsonWrapper writeToJSON() const;

  /**
   * @brief Read from the stream the serialized ciphertext (`Ctxt`) object using
   * JSON format.
   * @param str Input `std::istream`.
   * @param pubKey The `PubKey` to be used.
   * @return The deserialized `Ctxt` object.
   **/
  static Ctxt readFromJSON(std::istream& str, const PubKey& pubKey);

  /**
   * @brief Read from the `JsonWrapper` the serialized ciphertext (`Ctxt`)
   * object.
   * @param j The `JsonWrapper` containing the serialized `Ctxt` object.
   * @param pubKey The `PubKey` to be used.
   * @return The deserialized `Ctxt` object.
   **/
  static Ctxt readFromJSON(const JsonWrapper& j, const PubKey& pubKey);

  /**
   * @brief In-place read from the `str` `std::istream` the serialized
   * ciphertext (`Ctxt`) object.
   * @param j The `JsonWrapper` containing the serialized `Ctxt` object.
   **/
  void readJSON(std::istream& str);

  /**
   * @brief In-place read from the `JsonWrapper` the serialized ciphertext
   * (`Ctxt`) object.
   * @param j The `JsonWrapper` containing the serialized `Ctxt` object.
   **/
  void readJSON(const JsonWrapper& j);

  // scale up c1, c2 so they have the same ratFactor
  static void equalizeRationalFactors(Ctxt& c1, Ctxt& c2);

  // adds noise for CKKS deryption to mitigate against
  // plaintext leakage attacks, as in the paper
  // "On the Security of Homomorphic Encryption on Approximate Numbers",
  // by Li and Micciancio.
  void addedNoiseForCKKSDecryption(const SecKey& sk,
                                   double eps,
                                   NTL::ZZX& noise) const;
};

// set out=prod_{i=0}^{n-1} v[j], takes depth log n and n-1 products
// out could point to v[0], but having it pointing to any other v[i]
// will make the result unpredictable.
void totalProduct(Ctxt& out, const std::vector<Ctxt>& v);

//! For i=n-1...0, set v[i]=prod_{j<=i} v[j]
//! This implementation uses depth log n and (nlog n)/2 products
void incrementalProduct(std::vector<Ctxt>& v);

void innerProduct(Ctxt& result,
                  const std::vector<Ctxt>& v1,
                  const std::vector<Ctxt>& v2);
inline Ctxt innerProduct(const std::vector<Ctxt>& v1,
                         const std::vector<Ctxt>& v2)
{
  Ctxt ret(v1[0].getPubKey());
  innerProduct(ret, v1, v2);
  return ret;
}

//! Compute the inner product of a vectors of ciphertexts and a constant vector
void innerProduct(Ctxt& result,
                  const std::vector<Ctxt>& v1,
                  const std::vector<DoubleCRT>& v2);
inline Ctxt innerProduct(const std::vector<Ctxt>& v1,
                         const std::vector<DoubleCRT>& v2)
{
  Ctxt ret(v1[0].getPubKey());
  innerProduct(ret, v1, v2);
  return ret;
}

void innerProduct(Ctxt& result,
                  const std::vector<Ctxt>& v1,
                  const std::vector<NTL::ZZX>& v2);
inline Ctxt innerProduct(const std::vector<Ctxt>& v1,
                         const std::vector<NTL::ZZX>& v2)
{
  Ctxt ret(v1[0].getPubKey());
  innerProduct(ret, v1, v2);
  return ret;
}

//! frobeniusAutomorph: free function version of frobeniusAutomorph method
inline void frobeniusAutomorph(Ctxt& ctxt, long j)
{
  ctxt.frobeniusAutomorph(j);
}

//! conjugate: free function that is equivalent to frobeniusAutomorph(ctxt, 1)
inline void conjugate(Ctxt& ctxt) { frobeniusAutomorph(ctxt, 1); }

//! Free functions to extract real and imaginary parts.
//! Unlike the corresponding Ctxt method, extractImPart is thread safe.
//! NOTES: (1) these will raise an error for BGV ciphertexts.
void extractRealPart(Ctxt& c);
void extractImPart(Ctxt& c);

//! power: free function version of power method
inline void power(Ctxt& ctxt, long e) { ctxt.power(e); }

//! negate: free function version of negate method
inline void negate(Ctxt& ctxt) { ctxt.negate(); }

//! print to cerr some info about ciphertext
void CheckCtxt(const Ctxt& c, const char* label);

/**
 * @brief Extract the mod-p digits of a mod-p^r ciphertext.
 *
 * extractDigits returns in the slots of digits[j] the j'th-lowest digits
 * from the integers in the slots of the input. Namely, the i'th slot of
 * digits[j] contains the j'th digit in the p-base expansion of the integer
 * in the i'th slot of the *this.
 *
 * If r==0 then it is set to c.effectiveR(). It is assumed that the slots
 * of *this contains integers mod p^r, i.e., that only the free terms are
 * nonzero. If that assumptions does not hold then the result will not be
 * a valid ciphertext anymore.
 *
 * The "shortcut" flag is deprecated, it often leads to catastrophic failure
 * in the noise estimate. Calling the function with shortcut=true has not
 * effect, except printing a warning message to cerr.
 *
 * The output ciphertext digits[j] contains the j'th digit in the base-p
 * expansion of the input, and its plaintext space is modulo p^{r-j}. All
 * the ciphertexts in the output are at the same level.
 **/
void extractDigits(std::vector<Ctxt>& digits, const Ctxt& c, long r = 0);
// implemented in extractDigits.cpp

inline void extractDigits(std::vector<Ctxt>& digits,
                          const Ctxt& c,
                          long r,
                          bool shortCut)
{
  if (shortCut)
    std::cerr << "extractDigits: the shortCut flag is disabled\n";
  extractDigits(digits, c, r);
}

inline void Ctxt::extractBits(std::vector<Ctxt>& bits, long nBits2extract)
{
  extractDigits(bits, *this, nBits2extract);
}

// @brief Extract the mod-p digits of a mod-p^{r+e} ciphertext.

// extendExtractDigits assumes that the slots of *this contains integers mod
// p^{r+e} i.e., that only the free terms are nonzero. (If that assumptions
// does not hold then the result will not be a valid ciphertext anymore.)
//
// It returns in the slots of digits[j] the j'th-lowest digits from the
// integers in the slots of the input. Namely, the i'th slot of digits[j]
// contains the j'th digit in the p-base expansion of the integer in the i'th
// slot of the *this.  The plaintext space of digits[j] is mod p^{e+r-j}.

void extendExtractDigits(std::vector<Ctxt>& digits,
                         const Ctxt& c,
                         long r,
                         long e);
// implemented in extractDigits.cpp

} // namespace helib

#endif // ifndef HELIB_CTXT_H
