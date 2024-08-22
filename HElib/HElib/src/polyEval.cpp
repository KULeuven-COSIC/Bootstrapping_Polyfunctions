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
// Notice: this file was modified from HElib
#include <helib/Context.h>
#include <helib/polyEval.h>

#include <sstream>
#include <algorithm>
#include <math.h>
#include <cmath>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

namespace helib {

#include "Magma_declare.h"

// Returns the e'th power of X, computing it as needed
Ctxt& DynamicCtxtPowers::getPower(long e)
{
  if (v.at(e - 1).isEmpty()) { // Not computed yet, compute it now

    // largest power of two smaller than e
    long k = 1L << (NTL::NextPowerOfTwo(e) - 1);
    v[e - 1] = getPower(e - k); // compute X^e = X^{e-k} * X^k
    v[e - 1].multiplyBy(getPower(k));
    // FIXME: could drop down / cleanup further as an optimization?
  }
  return v[e - 1];
}

// Local functions for polynomial evaluation in some special cases
static void simplePolyEval(Ctxt& ret,
                           const NTL::ZZX& poly,
                           DynamicCtxtPowers& babyStep);
static void PatersonStockmeyer(Ctxt& ret,
                               const NTL::ZZX& poly,
                               long k,
                               long t,
                               long delta,
                               DynamicCtxtPowers& babyStep,
                               DynamicCtxtPowers& giantStep);
static void degPowerOfTwo(Ctxt& ret,
                          const NTL::ZZX& poly,
                          long k,
                          DynamicCtxtPowers& babyStep,
                          DynamicCtxtPowers& giantStep);
static void recursivePolyEval(Ctxt& ret,
                              const NTL::ZZX& poly,
                              long k,
                              DynamicCtxtPowers& babyStep,
                              DynamicCtxtPowers& giantStep);

static void recursivePolyEval(Ctxt& ret,
                              const Ctxt poly[],
                              long nCoeffs,
                              const NTL::Vec<Ctxt>& powers);

// Main entry point: Evaluate an encrypted polynomial on an encrypted input
// return in ret = sum_i poly[i] * x^i
void polyEval(Ctxt& ret, const NTL::Vec<Ctxt>& poly, const Ctxt& x)
{
  if (poly.length() <= 1) { // Some special cases
    if (poly.length() == 0)
      ret.clear(); // empty polynomial
    else
      ret = poly[0]; // constant polynomial
    return;
  }
  long deg = poly.length() - 1;

  long logD = NTL::NextPowerOfTwo(divc(poly.length(), 3));
  long d = 1L << logD;

  // We have d <= deg(poly) < 3d
  assertInRange(deg, d, 3l * d, "Poly degree not in [d, 3d)");

  NTL::Vec<Ctxt> powers(NTL::INIT_SIZE, logD + 1, x);
  if (logD > 0) {
    powers[1].square();
    for (long i = 2; i <= logD; i++) { // powers[i] = x^{2^i}
      powers[i] = powers[i - 1];
      powers[i].square();
    }
  }

  // Compute in three parts p0(X) + ( p1(X) + p2(X)*X^d )*X^d
  Ctxt tmp(ZeroCtxtLike, ret);
  recursivePolyEval(ret,
                    &poly[d],
                    std::min(d, poly.length() - d),
                    powers); // p1(X)

  if (poly.length() > 2 * d) { // p2 is not empty
    recursivePolyEval(tmp,
                      &poly[2 * d],
                      poly.length() - 2 * d,
                      powers); // p2(X)
    tmp.multiplyBy(powers[logD]);
    ret += tmp;
  }
  ret.multiplyBy(powers[logD]); // ( p1(X) + p2(X)*X^d )*X^d

  recursivePolyEval(tmp, &poly[0], d, powers); // p0(X)
  ret += tmp;
}

static void recursivePolyEval(Ctxt& ret,
                              const Ctxt poly[],
                              long nCoeffs,
                              const NTL::Vec<Ctxt>& powers)
{
  if (nCoeffs <= 1) { // edge condition
    if (nCoeffs == 0)
      ret.clear(); // empty polynomial
    else
      ret = poly[0]; // constant polynomial
    return;
  }
  long logD = NTL::NextPowerOfTwo(nCoeffs) - 1;
  long d = 1L << logD;
  Ctxt tmp(ZeroCtxtLike, ret);
  recursivePolyEval(tmp, &(poly[d]), nCoeffs - d, powers);
  recursivePolyEval(ret, &(poly[0]), d, powers);
  tmp.multiplyBy(powers[logD]);
  ret += tmp;
}

// Main entry point: Evaluate a cleartext polynomial on an encrypted input
void polyEval(Ctxt& ret, NTL::ZZX poly, const Ctxt& x, long k)
// Note: poly is passed by value, so caller keeps the original
{
  if (deg(poly) <= 2) {  // nothing to optimize here
    if (deg(poly) < 1) { // A constant
      ret.clear();
      ret.addConstant(coeff(poly, 0));
    } else { // A linear or quadratic polynomial
      DynamicCtxtPowers babyStep(x, deg(poly));
      simplePolyEval(ret, poly, babyStep);
    }
    return;
  }

  // How many baby steps: set k~sqrt(n/2), rounded up/down to a power of two

  // FIXME: There may be some room for optimization here: it may be possible
  // to choose k as something other than a power of two and still maintain
  // optimal depth, in principle we can try all possible values of k between
  // two consecutive powers of two and choose the one that gives the least
  // number of multiplies, conditioned on minimum depth.

  if (k <= 0) {
    long kk = (long)sqrt(deg(poly) / 2.0);
    k = 1L << NTL::NextPowerOfTwo(kk);

    // heuristic: if k>>kk then use a smaller power of two
    if ((k == 16 && deg(poly) > 167) || (k > 16 && k > (1.44 * kk)))
      k /= 2;
  }
#ifdef HELIB_DEBUG
  std::cerr << "  k=" << k;
#endif

  long n = divc(deg(poly), k); // n = ceil(deg(p)/k), deg(p) >= k*n
  DynamicCtxtPowers babyStep(x, k);
  const Ctxt& x2k = babyStep.getPower(k);

  // Special case when deg(p)>k*(2^e -1)
  if (n == (1L << NTL::NextPowerOfTwo(n))) { // n is a power of two
    DynamicCtxtPowers giantStep(x2k, n / 2);
    degPowerOfTwo(ret, poly, k, babyStep, giantStep);
    return;
  }

  // If n is not a power of two, ensure that poly is monic and that
  // its degree is divisible by k, then call the recursive procedure

  const NTL::ZZ p = NTL::to_ZZ(x.getPtxtSpace());
  NTL::ZZ top = LeadCoeff(poly);
  NTL::ZZ topInv; // the inverse mod p of the top coefficient of poly (if any)
  bool divisible = (n * k == deg(poly)); // is the degree divisible by k?
  long nonInvertible = InvModStatus(topInv, top, p);
  // 0 if invertible, 1 if not

  // FIXME: There may be some room for optimization below: instead of
  // adding a term X^{n*k} we can add X^{n'*k} for some n'>n, so long
  // as n' is smaller than the next power of two. We could save a few
  // multiplications since giantStep[n'] may be easier to compute than
  // giantStep[n] when n' has fewer 1's than n in its binary expansion.

  // extra!=0 denotes an added term extra*X^{n*k}
  NTL::ZZ extra = NTL::ZZ::zero();
  if (!divisible || nonInvertible) { // need to add a term
    top = NTL::to_ZZ(1);             // new top coefficient is one
    topInv = top;                    // also the new inverse is one
    // set extra = 1 - current-coeff-of-X^{n*k}
    extra = SubMod(top, coeff(poly, n * k), p);
    SetCoeff(poly, n * k); // set the top coefficient of X^{n*k} to one
  }

  long t = IsZero(extra) ? divc(n, 2) : n;
  DynamicCtxtPowers giantStep(x2k, t);

  if (!IsOne(top)) {
    poly *= topInv; // Multiply by topInv to make into a monic polynomial
    for (long i = 0; i <= n * k; i++)
      rem(poly[i], poly[i], p);
    poly.normalize();
  }
  recursivePolyEval(ret, poly, k, babyStep, giantStep);

  if (!IsOne(top)) {
    ret.multByConstant(top);
  }

  if (!IsZero(extra)) { // if we added a term, now is the time to subtract back
    Ctxt topTerm = giantStep.getPower(n);
    topTerm.multByConstant(extra);
    ret -= topTerm;
  }
}

// Simple evaluation sum f_i * X^i, assuming that babyStep has enough powers
static void simplePolyEval(Ctxt& ret,
                           const NTL::ZZX& poly,
                           DynamicCtxtPowers& babyStep)
{
  ret.clear();
  if (deg(poly) < 0)
    return; // the zero polynomial always returns zero

  // Ensure that we have enough powers

  // ensure that we have enough powers
  assertTrue(deg(poly) <= babyStep.size(),
             "BabyStep has not enough powers "
             "(required more than deg(poly))");

  NTL::ZZ coef;
  NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
  for (long i = 1; i <= deg(poly); i++) {
    rem(coef, coeff(poly, i), p);
    if (coef > p / 2)
      coef -= p;

    Ctxt tmp = babyStep.getPower(i); // X^i
    tmp.multByConstant(coef);        // f_i X^i
    ret += tmp;
  }
  // Add the free term
  rem(coef, ConstTerm(poly), p);
  if (coef > p / 2)
    coef -= p;
  ret.addConstant(coef);
  //  if (verbose) checkPolyEval(ret, babyStep[0], poly);
}

// The recursive procedure in the Paterson-Stockmeyer
// polynomial-evaluation algorithm from SIAM J. on Computing, 1973.
// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains >= k+delta powers
static void PatersonStockmeyer(Ctxt& ret,
                               const NTL::ZZX& poly,
                               long k,
                               long t,
                               long delta,
                               DynamicCtxtPowers& babyStep,
                               DynamicCtxtPowers& giantStep)
{
  if (deg(poly) <= babyStep.size()) { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }
  NTL::ZZX r = trunc(poly, k * t);      // degree <= k*2^e-1
  NTL::ZZX q = RightShift(poly, k * t); // degree == k(2^e-1) +delta

  const NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
  const NTL::ZZ& coef = coeff(r, deg(q));
  SetCoeff(r, deg(q), coef - 1); // r' = r - X^{deg(q)}

  NTL::ZZX c, s;
  DivRem(c, s, r, q); // r' = c*q + s
  // deg(s)<deg(q), and if c!= 0 then deg(c)<k-delta

  assertTrue(deg(s) < deg(q), "Degree of s is not less than degree of q");
  assertTrue(IsZero(c) || deg(c) < k - delta,
             "Nonzero c has not degree smaller than k - delta");
  SetCoeff(s, deg(q)); // s' = s + X^{deg(q)}, deg(s)==deg(q)

  // reduce the coefficients modulo p
  for (long i = 0; i <= deg(c); i++)
    rem(c[i], c[i], p);
  c.normalize();
  for (long i = 0; i <= deg(s); i++)
    rem(s[i], s[i], p);
  s.normalize();

  // Evaluate recursively poly = (c+X^{kt})*q + s'
  PatersonStockmeyer(ret, q, k, t / 2, delta, babyStep, giantStep);

  Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
  simplePolyEval(tmp, c, babyStep);
  tmp += giantStep.getPower(t);
  ret.multiplyBy(tmp);

  PatersonStockmeyer(tmp, s, k, t / 2, delta, babyStep, giantStep);
  ret += tmp;
}

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains >= k + (deg(poly) mod k) powers
static void degPowerOfTwo(Ctxt& ret,
                          const NTL::ZZX& poly,
                          long k,
                          DynamicCtxtPowers& babyStep,
                          DynamicCtxtPowers& giantStep)
{
  if (deg(poly) <= babyStep.size()) { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }
  long n = deg(poly) / k;                     // We assume n=2^e or n=2^e -1
  n = 1L << NTL::NextPowerOfTwo(n);           // round up to n=2^e
  NTL::ZZX r = trunc(poly, (n - 1) * k);      // degree <= k(2^e-1)-1
  NTL::ZZX q = RightShift(poly, (n - 1) * k); // 0 < degree < 2k
  SetCoeff(r, (n - 1) * k);                   // monic, degree == k(2^e-1)
  q -= 1;

  PatersonStockmeyer(ret, r, k, n / 2, 0, babyStep, giantStep);

  Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
  simplePolyEval(tmp, q, babyStep); // evaluate q

  // multiply by X^{k(n-1)} with minimum depth
  for (long i = 1; i < n; i *= 2) {
    tmp.multiplyBy(giantStep.getPower(i));
  }
  ret += tmp;
}

static void recursivePolyEval(Ctxt& ret,
                              const NTL::ZZX& poly,
                              long k,
                              DynamicCtxtPowers& babyStep,
                              DynamicCtxtPowers& giantStep)
{
  if (deg(poly) <= babyStep.size()) { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }

  long delta = deg(poly) % k;              // deg(poly) mod k
  long n = divc(deg(poly), k);             // ceil( deg(poly)/k )
  long t = 1L << (NTL::NextPowerOfTwo(n)); // t >= n, so t*k >= deg(poly)

  // Special case for deg(poly) = k * 2^e +delta
  if (n == t) {
    degPowerOfTwo(ret, poly, k, babyStep, giantStep);
    return;
  }

  // When deg(poly) = k*(2^e -1) we use the Paterson-Stockmeyer recursion
  if (n == t - 1 && delta == 0) {
    PatersonStockmeyer(ret, poly, k, t / 2, delta, babyStep, giantStep);
    return;
  }

  t = t / 2;

  // In any other case we have kt < deg(poly) < k(2t-1). We then set
  // u = deg(poly) - k*(t-1) and poly = q*X^u + r with deg(r)<u
  // and recurse on poly = (q-1)*X^u + (X^u+r)

  long u = deg(poly) - k * (t - 1);
  NTL::ZZX r = trunc(poly, u);      // degree <= u-1
  NTL::ZZX q = RightShift(poly, u); // degree == k*(t-1)
  q -= 1;
  SetCoeff(r, u); // degree == u

  PatersonStockmeyer(ret, q, k, t / 2, 0, babyStep, giantStep);

  Ctxt tmp = giantStep.getPower(u / k);
  if (delta != 0) { // if u is not divisible by k then compute it
    tmp.multiplyBy(babyStep.getPower(delta));
  }
  ret.multiplyBy(tmp);

  recursivePolyEval(tmp, r, k, babyStep, giantStep);
  ret += tmp;
}

// raise ciphertext to some power
void Ctxt::power(long e)
{
  if (e < 1)
    throw InvalidArgument("Cannot raise a ctxt to a non positive exponent");

  if (e == 1)
    return; // nothing to do

  long ell = NTL::NumBits(e); // e < 2^l <= 2e

  if (static_cast<unsigned long>(e) ==
      (1UL << (ell - 1))) { // e is a power of two, just square enough times
    while (--ell > 0)
      square();
    return;
  }

  // Otherwise use the "DynamicCtxtPowers" from polyEval, it uses e Ctxt
  // objects as temporary space but keeps the level as low as possible
  DynamicCtxtPowers pwrs(*this, e);
  *this = pwrs.getPower(e);
}

long ceiling(double number) {
    return lround(number + 0.49999);
}

long flooring(double number) {
    return lround(number - 0.49999);
}

long floorPowerOfTwo(long number) {
    return std::pow(2, flooring(log(number) / log(2)));
}

bool isPowerOfTwo(long number) {
    return floorPowerOfTwo(number) == number;
}

// Compute the common spacing of the given polynomials, i.e., the largest value of n
// such that we have poly = f(x^n) for some polynomial f(x)
long getSpacing(const std::vector<NTL::ZZX>& poly_list) {
    long spacing(deg(poly_list[0]));
    for (NTL::ZZX poly : poly_list) {
        long last_non_zero_index = 0;
        for (long index = 1; index <= deg(poly); index++) {
            if (coeff(poly, index) != 0) {
                spacing = NTL::GCD(spacing, index - last_non_zero_index);
                last_non_zero_index = index;
            }
        }
    }
    return spacing;
}

// Check if the given sequence of polynomials is odd
bool areOddPolynomials(const std::vector<NTL::ZZX>& poly_list) {
    // Check if each individual polynomial is odd
    for (NTL::ZZX poly : poly_list)
        for (long index = 0; index <= deg(poly); index++)
            if ((index % 2 == 0) && (coeff(poly, index) != 0))
                return false;

    return true;
}

// Return the parameters that lead to the smallest number of non-constant multiplications
// Degree of the polynomials is at least k * (2 ^ m)
PS_parameters getBestParameters(const std::vector<NTL::ZZX>& polynomials, bool lazy) {
    for (NTL::ZZX polynomial : polynomials) {
        assertNeq(deg(polynomial), (long)-1, "Degree should be positive.");
        assertNeq(deg(polynomial), (long)0, "Degree should be positive.");
    }

    // Maximum degree of the given polynomials
    long d = 0;
    for (NTL::ZZX polynomial : polynomials)
        d = std::max(d, deg(polynomial));

    // Check if all polynomials are odd
    bool odd = areOddPolynomials(polynomials);

    // Compute best set of parameters by iteration
    int bestM = 0;
    int bestK = 0;
    int bestMultiplications = -1;
    bool bestOdd = false;
    for (int m = 0; m <= ceiling(log(d) / log(2)); m++) {
        // Compute corresponding k parameter and number of multiplications (start with baby step only)
        // Note that we cannot combine lazy rescaling with odd polynomials (different computation in the baby step)
        // --> Lazy rescaling is prioritized since it can be set as a flag in the parameter list
        // Edit: the above is not true, but it was not updated in the implementation for reproducibility with the paper
        int k = ceiling(d / pow(2, m));
        int nbMultiplications;
        bool currentOdd = false;
        if (lazy) {
            if (k == 1)
                nbMultiplications = m - 1;
            else
                nbMultiplications = ((k - 1) / 2) + m;
        } else {
            if (m == 0)
                nbMultiplications = k - 1;
            else
                nbMultiplications = k + m - 2;

            // Possibly use different algorithm if polynomials are odd (only if operation count is better)
            if (odd) {
                int k_odd = k;
                int nbMultiplicationsOdd;
                if (m == 0) {
                    nbMultiplicationsOdd = (k_odd / 2) + flooring(log(k_odd) / log(2));
                } else {
                    // The optimized procedure only works for even k
                    if ((k_odd % 2) == 1)
                        k_odd += 1;

                    // Make sure that we can always compute x^k as a product of two factors that were computed in the baby step
                    // This is done by multiplying x^e and x^d where either both e and d are odd or both are a power of 2
                    // Note that this is not always possible if k == 0 (mod 4) and we might have to increase k by 2
                    // Edit: the above is not true, but it was not updated in the implementation for reproducibility with the paper
                    int remaining_exponent = k_odd - floorPowerOfTwo(k_odd - 1);
                    if (((k_odd % 4) == 0) && (!isPowerOfTwo(remaining_exponent)))
                        k_odd += 2;

                    // Now compute the actual number of multiplications
                    nbMultiplicationsOdd = (k_odd / 2) + flooring(log(k_odd - 1) / log(2)) + m - 1;
                }

                // Check if odd parameters are better
                if (nbMultiplicationsOdd < nbMultiplications) {
                    k = k_odd;
                    nbMultiplications = nbMultiplicationsOdd;
                    currentOdd = true;
                }
            }
        }

        // Add extra number for giant step
        for (NTL::ZZX polynomial : polynomials) {
            nbMultiplications += (ceiling(deg(polynomial) / k) - 1);
            if (lazy) {    // One extra non-scalar multiplication in giant step
                nbMultiplications += 1;
                int deg_mod = deg(polynomial) % k;    // One less non-scalar multiplication if baby step has only linear terms
                if ((deg_mod != 0) && (deg_mod <= ((k + 1) / 2)))
                    nbMultiplications -= 1;
            }
        }

        // Check whether the parameters are better than the current best ones
        if ((bestMultiplications == -1) || (nbMultiplications < bestMultiplications)) {
            bestM = m;
            bestK = k;
            bestMultiplications = nbMultiplications;
            bestOdd = currentOdd;
        }
    }
    PS_parameters result;
    result.m = bestM; result.k = bestK; result.multiplications = bestMultiplications; result.odd = bestOdd;
    return result;
}

// Preprocessing for polynomial evaluation
// Evaluate x^spacing and compute updated polynomials accordingly
void customPolyEvalPreprocessing(Ctxt& element, const std::vector<NTL::ZZX>& input, std::vector<NTL::ZZX>& output) {
    long spacing = getSpacing(input);
    element.power(spacing);
    for (NTL::ZZX poly : input) {
        output.push_back(NTL::ZZX());
        for (long index = 0; index <= deg(poly) / spacing; index++)
            SetCoeff(output.back(), index, coeff(poly, index * spacing));
    }
}

// Recursive part of Paterson-Stockmeyer algorithm
void customPolyEvalRecursive(Ctxt& result, const std::vector<NTL::ZZ>& coeff, std::vector<Ctxt>& xExp1, std::vector<Ctxt>& xExp2, int m, int k, bool lazy) {
    // Base cases
    result = Ctxt(ZeroCtxtLike, xExp1[0]);
    if (coeff.size() == 0) {
        return;
    } else if (m == 0) {
        for (std::vector<NTL::ZZ>::size_type index = 0; index < coeff.size(); index++) {  // Inner loop: baby step
            if (coeff[index] != 0) {
                Ctxt tmp(xExp1[index]);
                if (coeff[index] == -1)
                    tmp.negate();
                else if (coeff[index] != 1)
                    tmp.multByConstant(coeff[index]);
                result.addCtxt(tmp);
            }
        }
        return;
    }

    // Recursive case
    Ctxt tmp(ZeroCtxtLike, xExp1[0]);
    long index = std::min<long>(k * pow(2, m - 1), coeff.size());
    customPolyEvalRecursive(result, std::vector<NTL::ZZ>(coeff.begin(), coeff.begin() + index), xExp1, xExp2, m - 1, k, lazy);
    customPolyEvalRecursive(tmp, std::vector<NTL::ZZ>(coeff.begin() + index, coeff.end()), xExp1, xExp2, m - 1, k, lazy);
    tmp.customMultiplyBy(xExp2[m - 1], lazy);
    result.addCtxt(tmp);
}

// Evaluate the given polynomials in the given element: the algorithm is optimized for lowest number of
// multiplications since the depth is already optimal (counting only non-scalar multiplications)
// This function can also execute the lazy baby-step/giant-step algorithm if the lazy flag is set to true
void customPolyEval(std::vector<Ctxt>& result, const std::vector<NTL::ZZX>& polynomials, const Ctxt& element, bool lazy) {
    for (NTL::ZZX polynomial : polynomials) {
        assertNeq(deg(polynomial), (long)-1, "Degree should be positive.");
        assertNeq(deg(polynomial), (long)0, "Degree should be positive.");
    }

    // Evaluate x^spacing and update polynomials accordingly
    // Also determine the optimal parameters for the remaining polynomials
    Ctxt new_element(element);
    std::vector<NTL::ZZX> new_polynomials;
    customPolyEvalPreprocessing(new_element, polynomials, new_polynomials);
    PS_parameters parameters = getBestParameters(new_polynomials, lazy);

    // Precompute x ^ exp with exp = 1, ..., k
    std::vector<Ctxt> xExp1{new_element};
    for (int exp = 2; exp <= parameters.k; exp++) {
        int ind1, ind2;
        if (parameters.odd) {   // For odd polynomials, we use the algorithm that only adapts the baby step (not the one that rewrites to x*f(x^2) because of depth increase)
            if ((exp % 2) == 0) {
                if (isPowerOfTwo(exp) || (exp == parameters.k)) {
                    ind1 = (((exp % 4) == 0) ? floorPowerOfTwo(exp - 1) : (exp / 2));
                    ind2 = exp - ind1;
                    Ctxt tmp(xExp1[ind1 - 1]);
                    tmp.customMultiplyBy(xExp1[ind2 - 1], lazy);
                    xExp1.push_back(tmp);
                } else {
                    xExp1.push_back(Ctxt(ZeroCtxtLike, element));   // Just append garbage
                }
            } else {
                ind1 = floorPowerOfTwo(exp);
                ind2 = exp - ind1;
                Ctxt tmp(xExp1[ind1 - 1]);
                tmp.customMultiplyBy(xExp1[ind2 - 1], lazy);
                xExp1.push_back(tmp);
            }
        } else {
            // Choose indices such that the depth is as low as possible
            ind1 = exp / 2;
            ind2 = exp - ind1;

            xExp1[ind1 - 1].reLinearize();
            xExp1[ind2 - 1].reLinearize();
            Ctxt tmp(xExp1[ind1 - 1]);
            tmp.customMultiplyBy(xExp1[ind2 - 1], lazy);
            xExp1.push_back(tmp);
        }
    }

    // Sanitize result for giant step
    if (parameters.m != 0)
        xExp1.back().reLinearize();

    // Precompute x ^ exp with exp = k, 2 * k, ..., (2 ^ (m - 1)) * k
    std::vector<Ctxt> xExp2{xExp1.back()};
    for (int exp = 1; exp < parameters.m; exp++) {
        Ctxt tmp(xExp2.back());
        tmp.multiplyBy(tmp);
        xExp2.push_back(tmp);
    }

    // Compute evaluation for each of the polynomials
    result.clear();
    for (NTL::ZZX polynomial : new_polynomials) {
        // Compute coefficients as a list
        std::vector<NTL::ZZ> coeff_list;
        for (int exp = 0; exp <= deg(polynomial); exp++) {
            coeff_list.push_back(coeff(polynomial, exp));
        }

        // Return result via sequence of recursive calls
        result.push_back(Ctxt(ZeroCtxtLike, new_element));
        customPolyEvalRecursive(result.back(), std::vector<NTL::ZZ>(coeff_list.begin() + 1, coeff_list.end()), xExp1, xExp2, parameters.m, parameters.k, lazy);
        if (coeff_list[0] != 0)
            result.back().addConstant(coeff_list[0]);
        result.back().reLinearize();
    }
}

#if 0
/**********************************************************************/
/*     FOR DEBUGGING PURPOSES, the same procedure for plaintext x     */
/**********************************************************************/

static long verbose = 0;
class DynamicPtxtPowers {
private:
  long p;  // the modulus
  vec_long v;    // A vector storing the powers themselves
  vec_long dpth; // keeping track of depth for debugging purposes

public:
  DynamicPtxtPowers(long _x, long _p, long nPowers, long _d=1) : p(_p)
  {
    // Sanity check
    assertTrue<InvalidArgument>(_x >= 0l, "_x must be greater equal than 0");
    assertTrue<InvalidArgument>(_p > 1l, "_p must be greater than 1"); // Sanity check
    assertTrue<InvalidArgument>(nPowers > 0l, "nPowers must be greater than 0"); // Sanity check
    v.SetLength(nPowers);
    dpth.SetLength(nPowers);
    for (long i=1; i<nPowers; i++) // Initializes nPowers empty slots
      dpth[i]=v[i]=-1;

    v[0] = _x % p;         // store X itself in v[0]
    dpth[0] = _d;
    //    std::cerr << " initial power="<<v[0]<<", initial depth="<<dpth[0];
  }

  // Returns the e'th power, computing it as needed
  long getPower(long e); // must use e >= 1, else throws an exception

  // dp.at(i) and dp[i] both return the i+1st power
  long at(long i) { return getPower(i+1); }
  long operator[](long i) { return getPower(i+1); }

  const vec_long& getVector() const { return v; }
  long size() const { return v.length(); }
  bool wasComputed(long i) { return (i>=0 && i<v.length() && v[i]>=0); }
  long getDepth(long e) { return (e>0)? dpth.at(e-1) : 0; }
};

// Returns the e'th power, computing it as needed
long DynamicPtxtPowers::getPower(long e)
{
  // FIXME: Do we want to allow the vector to grow? If so then begin by
  // checking e<v.length() and resizing if not. Currently throws an exception.

  if (v.at(e-1)<0) { // Not computed yet, compute it now

    long k = 1L<<(NTL::NextPowerOfTwo(e)-1); // largest power of two smaller than e

    v[e-1] = getPower(e-k);             // compute X^e = X^{e-k} * X^k
    v[e-1] = MulMod(v[e-1], getPower(k), p);
    dpth[e-1] = max(getDepth(k),getDepth(e-k)) +1;
    nMults++;
  }
  return v[e-1];
}

std::ostream& operator<< (std::ostream &s, const DynamicPtxtPowers &d)
{
  return s << d.getVector();
}

static long
recursivePolyEval(const NTL::ZZX& poly, long k, DynamicPtxtPowers& babyStep,
		  DynamicPtxtPowers& giantStep, long mod,
		  long& recursiveDepth);
static long
PatersonStockmeyer(const NTL::ZZX& poly, long k, long t, long delta,
		   DynamicPtxtPowers& babyStep, DynamicPtxtPowers& giantStep,
		   long mod, long& recursiveDepth);

long simplePolyEval(const NTL::ZZX& poly, DynamicPtxtPowers& babyStep, long mod)
{
  // ensure that we have enough powers
  assertTrue(deg(poly)<=(long)babyStep.size(), "BabyStep has not enough powers (required more than deg(poly))");

  long ret = rem(ConstTerm(poly), mod);
  for (long i=0; i<deg(poly); i++) {
    long coeff = rem(poly[i+1], mod); // f_{i+1}
    long tmp = babyStep[i];           // X^{i+1}
    tmp = MulMod(tmp, coeff, mod);    // f_{i+1} X^{i+1}
    ret = AddMod(ret, tmp, mod);
  }
  return ret;
}

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains k+ (deg(poly) mod k) powers
static long degPowerOfTwo(const NTL::ZZX& poly, long k, DynamicPtxtPowers& babyStep,
			  DynamicPtxtPowers& giantStep, long mod,
			  long& recursiveDepth)
{
  if (deg(poly)<=babyStep.size()) { // Edge condition, use simple eval
    long ret = simplePolyEval(poly, babyStep, mod);
    recursiveDepth = babyStep.getDepth(deg(poly));
    return ret;
  }
  long subDepth1 =0, subDepth2=0;
  long n = deg(poly)/k;        // We assume n=2^e or n=2^e -1
  n = 1L << NTL::NextPowerOfTwo(n); // round up to n=2^e
  NTL::ZZX r = trunc(poly, (n-1)*k);      // degree <= k(2^e-1)-1
  NTL::ZZX q = RightShift(poly, (n-1)*k); // 0 < degree < 2k
  SetCoeff(r, (n-1)*k);              // monic, degree == k(2^e-1)
  q -= 1;
  if (verbose) std::cerr << ", recursing on "<<r<<" + X^"<<(n-1)*k<<"*"<<q<<std::endl;

  long ret = PatersonStockmeyer(r, k, n/2, 0,
				babyStep, giantStep, mod, subDepth2);
  if (verbose)
    std::cerr << "  PatersonStockmeyer("<<r<<") returns "<<ret
	 << ", depth="<<subDepth2<<std::endl;

  long tmp = simplePolyEval(q, babyStep, mod); // evaluate q
  subDepth1 = babyStep.getDepth(deg(q));
  if (verbose)
    std::cerr << "  simplePolyEval("<<q<<") returns "<<tmp
	 << ", depth="<<subDepth1<<std::endl;

  // multiply by X^{k(n-1)} with minimum depth
  for (long i=1; i<n; i*=2) {
    tmp = MulMod(tmp, giantStep.getPower(i), mod);
    nMults++;
    subDepth1 = max(subDepth1, giantStep.getDepth(i)) +1;
    if (verbose)
      std::cerr << "    after mult by giantStep.getPower("<<i<< ")="
	   << giantStep.getPower(i)<<" of depth="<< giantStep.getDepth(i)
	   << ",  ret="<<tmp<<" and depth is "<<subDepth1<<std::endl;
  }
  totalDepth = max(subDepth1, subDepth2);
  return AddMod(ret, tmp, mod); // return q * X^{k(n-1)} + r
}

// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains k+delta powers
static long
PatersonStockmeyer(const NTL::ZZX& poly, long k, long t, long delta,
		   DynamicPtxtPowers& babyStep, DynamicPtxtPowers& giantStep,
		   long mod, long& recursiveDepth)
{
  if (verbose) std::cerr << "PatersonStockmeyer("<<poly<<"), k="<<k<<", t="<<t<<std::endl;

  if (deg(poly)<=babyStep.size()) { // Edge condition, use simple eval
    long ret = simplePolyEval(poly, babyStep, mod);
    recursiveDepth = babyStep.getDepth(deg(poly));
    return ret;
  }
  long subDepth1=0, subDepth2=0;
  long ret, tmp;

  NTL::ZZX r = trunc(poly, k*t);      // degree <= k*2^e-1
  NTL::ZZX q = RightShift(poly, k*t); // degree == k(2^e-1) +delta

  if (verbose) std::cerr << "  r ="<<r<< ", q="<<q;

  const NTL::ZZ& coef = coeff(r,deg(q));
  SetCoeff(r, deg(q), coef-1);  // r' = r - X^{deg(q)}

  if (verbose) std::cerr << ", r'="<<r;

  NTL::ZZX c,s;
  DivRem(c,s,r,q); // r' = c*q + s
  // deg(s)<deg(q), and if c!= 0 then deg(c)<k-delta

  if (verbose) std::cerr << ", c="<<c<< ", s ="<<s<<std::endl;

  assertTrue(deg(s)<deg(q), "Degree of s is not less than degree of q");
  assertTrue(IsZero(c) || deg(c)<k - delta, "Nonzero c has not degree smaller than k - delta");

  SetCoeff(s,deg(q)); // s' = s + X^{deg(q)}, deg(s)==deg(q)

  // reduce the coefficients modulo mod
  for (long i=0; i<=deg(c); i++) rem(c[i],c[i], NTL::to_ZZ(mod));
  c.normalize();
  for (long i=0; i<=deg(s); i++) rem(s[i],s[i], NTL::to_ZZ(mod));
  s.normalize();

  if (verbose) std::cerr << " {t==n+1} recursing on "<<s<<" + (X^"<<t*k<<"+"<<c<<")*"<<q<<std::endl;

  // Evaluate recursively poly = (c+X^{kt})*q + s'
  tmp = simplePolyEval(c, babyStep, mod);
  tmp = AddMod(tmp, giantStep.getPower(t), mod);
  subDepth1 = max(babyStep.getDepth(deg(c)), giantStep.getDepth(t));

  ret = PatersonStockmeyer(q, k, t/2, delta,
			   babyStep, giantStep, mod, subDepth2);

  if (verbose) {
    std::cerr << "  PatersonStockmeyer("<<q<<") returns "<<ret
	 << ", depth="<<subDepth2<<std::endl;
    if (ret != polyEvalMod(q,babyStep[0], mod)) {
      std::stringstream ss;
      ss << "  **1st recursive call failed, q="<<q;
      throw RuntimeError(ss.get());
    }
  }
  ret = MulMod(ret, tmp, mod);
  nMults++;
  subDepth1 = max(subDepth1, subDepth2)+1;

  tmp = PatersonStockmeyer(s, k, t/2, delta,
			   babyStep, giantStep, mod, subDepth2);
  if (verbose) {
    std::cerr << "  PatersonStockmeyer("<<s<<") returns "<<tmp
	 << ", depth="<<subDepth2<<std::endl;
    if (tmp != polyEvalMod(s,babyStep[0], mod)) {
      std::stringstream ss;
      ss << "  **2nd recursive call failed, s="<<s;
      throw RuntimeError(ss.get());
    }
  }
  ret = AddMod(ret,tmp,mod);
  recursiveDepth = max(subDepth1, subDepth2);
  return ret;
}


// This procedure assumes that poly is monic and that babyStep contains
// at least k+delta powers, where delta = deg(poly) mod k
static long
recursivePolyEval(const NTL::ZZX& poly, long k, DynamicPtxtPowers& babyStep,
		  DynamicPtxtPowers& giantStep, long mod,
		  long& recursiveDepth)
{
  if (deg(poly)<=babyStep.size()) { // Edge condition, use simple eval
    long ret = simplePolyEval(poly, babyStep, mod);
    recursiveDepth = babyStep.getDepth(deg(poly));
    return ret;
  }

  if (verbose) std::cerr << "recursivePolyEval("<<poly<<")\n";

  long delta = deg(poly) % k; // deg(poly) mod k
  long n = divc(deg(poly),k); // ceil( deg(poly)/k )
  long t = 1L<<(NTL::NextPowerOfTwo(n)); // t >= n, so t*k >= deg(poly)

  // Special case for deg(poly) = k * 2^e +delta
  if (n==t)
    return degPowerOfTwo(poly, k, babyStep, giantStep, mod, recursiveDepth);

  // When deg(poly) = k*(2^e -1) we use the Paterson-Stockmeyer recursion
  if (n == t-1 && delta==0)
    return PatersonStockmeyer(poly, k, t/2, delta,
			      babyStep, giantStep, mod, recursiveDepth);

  t = t/2;

  // In any other case we have kt < deg(poly) < k(2t-1). We then set
  // u = deg(poly) - k*(t-1) and poly = q*X^u + r with deg(r)<u
  // and recurse on poly = (q-1)*X^u + (X^u+r)

  long u = deg(poly) - k*(t-1);
  NTL::ZZX r = trunc(poly, u);      // degree <= u-1
  NTL::ZZX q = RightShift(poly, u); // degree == k*(t-1)
  q -= 1;
  SetCoeff(r, u);              // degree == u

  long ret, tmp;
  long subDepth1=0, subDepth2=0;
  if (verbose)
    std::cerr << " {deg(poly)="<<deg(poly)<<"<k*(2t-1)="<<k*(2*t-1)
	 << "} recursing on "<<r<<" + X^"<<u<<"*"<<q<<std::endl;
  ret = PatersonStockmeyer(q, k, t/2, 0,
			   babyStep, giantStep, mod, subDepth1);

  if (verbose) {
    std::cerr << "  PatersonStockmeyer("<<q<<") returns "<<ret<<", depth="<<subDepth1<<std::endl;
    if (ret != polyEvalMod(q,babyStep[0], mod)) {
      std::stringstream ss;
      ss << "  @@1st recursive call failed, q="<<q
     	   << ", ret="<<ret<<"!=" << polyEvalMod(q,babyStep[0], mod);
      throw RuntimeError(ss.get());
    }
  }

  tmp = giantStep.getPower(u/k);
  subDepth2 = giantStep.getDepth(u/k);
  if (delta!=0) { // if u is not divisible by k then compute it
    if (verbose)
      std::cerr <<"  multiplying by X^"<<u
	   <<"=giantStep.getPower("<<(u/k)<<")*babyStep.getPower("<<delta<<")="
	   << giantStep.getPower(u/k)<<"*"<<babyStep.getPower(delta)
	   << "="<<tmp<<std::endl;
    tmp = MulMod(tmp, babyStep.getPower(delta), mod);
    nMults++;
    subDepth2++;
  }
  ret = MulMod(ret, tmp, mod);
  nMults ++;
  subDepth1 = max(subDepth1, subDepth2)+1;

  if (verbose) std::cerr << "  after mult by X^{k*"<<u<<"+"<<delta<<"}, depth="<< subDepth1<<std::endl;

  tmp = recursivePolyEval(r, k, babyStep, giantStep, mod, subDepth2);
  if (verbose)
    std::cerr << "  recursivePolyEval("<<r<<") returns "<<tmp<<", depth="<<subDepth2<<std::endl;
  if (tmp != polyEvalMod(r,babyStep[0], mod)) {
    std::stringstream ss;
    ss << "  @@2nd recursive call failed, r="<<r
      << ", ret="<<tmp<<"!=" << polyEvalMod(r,babyStep[0], mod);
    throw RuntimeError(ss.get());
  }
  recursiveDepth = max(subDepth1, subDepth2);
  return AddMod(ret, tmp, mod);
}

// Note: poly is passed by value, not by reference, so the calling routine
// keeps its original polynomial
long evalPolyTopLevel(NTL::ZZX poly, long x, long p, long k=0)
{
  if (verbose)
  std::cerr << "\n* evalPolyTopLevel: p="<<p<<", x="<<x<<", poly="<<poly;

  if (deg(poly)<=2) { // nothing to optimize here
    if (deg(poly)<1) return to_long(coeff(poly, 0));
    DynamicPtxtPowers babyStep(x, p, deg(poly));
    long ret = simplePolyEval(poly, babyStep, p);
    totalDepth = babyStep.getDepth(deg(poly));
    return ret;
  }

  // How many baby steps: set k~sqrt(n/2), rounded up/down to a power of two

  // FIXME: There may be some room for optimization here: it may be possible
  // to choose k as something other than a power of two and still maintain
  // optimal depth, in principle we can try all possible values of k between
  // the two powers of two and choose the one that goves the least number
  // of multiplies, conditioned on minimum depth.

  if (k<=0) {
    long kk = (long) sqrt(deg(poly)/2.0);
    k = 1L << NTL::NextPowerOfTwo(kk);

    // heuristic: if k>>kk then use a smaler power of two
    if ((k==16 && deg(poly)>167) || (k>16 && k>(1.44*kk)))
      k /= 2;
  }
  std::cerr << ", k="<<k;

  long n = divc(deg(poly),k);          // deg(p) = k*n +delta
  if (verbose) std::cerr << ", n="<<n<<std::endl;

  DynamicPtxtPowers babyStep(x, p, k);
  long x2k = babyStep.getPower(k);

  // Special case when deg(p)>k*(2^e -1)
  if (n==(1L << NTL::NextPowerOfTwo(n))) { // n is a power of two
    DynamicPtxtPowers giantStep(x2k, p, n/2, babyStep.getDepth(k));
    if (verbose)
      std::cerr << "babyStep="<<babyStep<<", giantStep="<<giantStep<<std::endl;
    long ret = degPowerOfTwo(poly, k, babyStep, giantStep, p, totalDepth);

    if (verbose) {
      std::cerr << "  degPowerOfTwo("<<poly<<") returns "<<ret<<", depth="<<totalDepth<<std::endl;
      if (ret != polyEvalMod(poly,babyStep[0], p)) {
        std::stringstream ss;
        ss << "  ## recursive call failed, ret="<<ret<<"!="
          << polyEvalMod(poly,babyStep[0], p);
        throw RuntimeError(ss.get());
      }
      // std::cerr << "  babyStep depth=[";
      // for (long i=0; i<babyStep.size(); i++)
      // 	std::cerr << babyStep.getDepth(i+1)<<" ";
      // std::cerr << "]\n";
      // std::cerr << "  giantStep depth=[";
      // for (long i=0; i<giantStep.size(); i++)
      // 	std::cerr<<giantStep.getDepth(i+1)<<" ";
      // std::cerr << "]\n";
    }
    return ret;
  }

  // If n is not a power of two, ensure that poly is monic and that
  // its degree is divisible by k, then call the recursive procedure

  NTL::ZZ topInv; // the inverse mod p of the top coefficient of poly (if any)
  bool divisible = (n*k == deg(poly)); // is the degree divisible by k?
  long nonInvertibe = InvModStatus(topInv, LeadCoeff(poly), NTL::to_ZZ(p));
       // 0 if invertible, 1 if not

  // FIXME: There may be some room for optimization below: instead of
  // adding a term X^{n*k} we can add X^{n'*k} for some n'>n, so long
  // as n' is smaller than the next power of two. We could save a few
  // multiplications since giantStep[n'] may be easier to compute than
  // giantStep[n] when n' has fewer 1's than n in its binary expansion.

  long extra = 0;        // extra!=0 denotes an added term extra*X^{n*k}
  if (!divisible || nonInvertibe) {  // need to add a term
    // set extra = 1 - current-coeff-of-X^{n*k}
    extra = SubMod(1, to_long(coeff(poly,n*k)), p);
    SetCoeff(poly, n*k); // set the top coefficient of X^{n*k} to one
    topInv = NTL::to_ZZ(1);   // inverse of new top coefficient is one
  }

  long t = (extra==0)? divc(n,2) : n;
  DynamicPtxtPowers giantStep(x2k, p, t, babyStep.getDepth(k));

  if (verbose)
    std::cerr << "babyStep="<<babyStep<<", giantStep="<<giantStep<<std::endl;

  long y; // the value to return
  long subDepth1 =0;
  if (!IsOne(topInv)) {
    long top = to_long(poly[n*k]); // record the current top coefficient
    //    std::cerr << ", top-coeff="<<top;

    // Multiply by topInv modulo p to make into a monic polynomial
    poly *= topInv;
    for (long i=0; i<=n*k; i++) rem(poly[i], poly[i], NTL::to_ZZ(p));
    poly.normalize();

    y = recursivePolyEval(poly, k, babyStep, giantStep, p, subDepth1);
    if (verbose) {
      std::cerr << "  recursivePolyEval("<<poly<<") returns "<<y<<", depth="<<subDepth1<<std::endl;
      if (y != polyEvalMod(poly,babyStep[0], p)) {
        std::stringstream ss;
        ss << "## recursive call failed, ret="<<y<<"!="
          << polyEvalMod(poly,babyStep[0], p);
        throw RuntimeError(ss.get());
      }
    }
    y = MulMod(y, top, p); // multiply by the original top coefficient
  }
  else {
    y = recursivePolyEval(poly, k, babyStep, giantStep, p, subDepth1);
    if (verbose) {
      std::cerr << "  recursivePolyEval("<<poly<<") returns "<<y<<", depth="<<subDepth1<<std::endl;
      if (y != polyEvalMod(poly,babyStep[0], p)) {
        std::stringstream ss;
        ss << "## recursive call failed, ret="<<y<<"!="
          << polyEvalMod(poly,babyStep[0], p);
        throw RuntimeError(ss.get());
      }
    }
  }

  if (extra != 0) { // if we added a term, now is the time to subtract back
    if (verbose) std::cerr << ", subtracting "<<extra<<"*X^"<<k*n;
    extra = MulMod(extra, giantStep.getPower(n), p);
    totalDepth = max(subDepth1, giantStep.getDepth(n));
    y = SubMod(y, extra, p);
  }
  else totalDepth = subDepth1;
  if (verbose) std::cerr << std::endl;
  return y;
}
#endif

// Read polynomial from file and store it in the given variable
void parsePolynomial(NTL::ZZX& poly_result, const char* fileName) {
  std::ifstream indata(fileName);
  assertTrue((bool)indata, "Polynomial file does not exist.");

  poly_result = NTL::ZZX();
  int index = 0;
  while (indata) {
    std::string tmp;
    indata >> tmp;
    if (!tmp.empty())
      SetCoeff(poly_result, index++, NTL::to_ZZ(tmp.c_str()));
  }
}

// Note that we store the polynomials in a global variable, which is bad practice
int init_polynomials() {
  for (long index = 0; index < nb_primes; index++) {
    for (long inner_index = 0; inner_index < nb_e_inner; inner_index++) {
      for (long e = e_inner_list[inner_index] + 1;; e++) {  // Load polynomials starting from precision = e_inner + 1
        std::stringstream stream1, stream2, stream3;
        stream1 << primes_list[index]; stream2 << e_inner_list[inner_index]; stream3 << e;
        std::string str1(stream1.str()), str2(stream2.str()), str3(stream3.str());

        std::string fileName((std::string("polynomials/poly") + str1 + std::string("_") + str2 + std::string("_") + str3 + std::string(".txt")).c_str());
        std::ifstream indata(fileName.c_str());
        if (!indata)
          break;

        NTL::ZZX result;
        parsePolynomial(result, fileName.c_str());
        polynomial_vector[index][inner_index].push_back(result);
      }
    }
  }
  return 0;
}

int nb_relin = 0;

#include "Magma_define.h"

static int ignore_this_var = init_polynomials();

} // namespace helib
