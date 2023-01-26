#include <iostream>

#include <helib/helib.h>
#include <helib/debugging.h>
#include <helib/polyEval.h>

#include <thread>
#include <chrono>

// Print the plaintext corresponding to the given ciphertext
void printPlaintext(const helib::Ctxt& ctxt, const helib::SecKey& secret_key, const char* text = "result") {
    long nslots = ctxt.getContext().getEA().size();
    std::vector<long> ptxt_res(nslots);
    ctxt.getContext().getEA().decrypt(ctxt, secret_key, ptxt_res);
    std::cout << "Decrypted " << text << ":" << std::endl;
    for (int i = 0; i < nslots; i++)
        std::cout << ptxt_res[i] << " ";
    std::cout << std::endl;
}

// Check if the given ciphertexts encrypt the same plaintext
bool areEqualCiphertexts(const helib::Ctxt& ctxt1, const helib::Ctxt& ctxt2, const helib::SecKey& secret_key) {
    long nslots = ctxt1.getContext().getEA().size();
    std::vector<long> ptxt_res1(nslots), ptxt_res2(nslots);
    ctxt1.getContext().getEA().decrypt(ctxt1, secret_key, ptxt_res1);
    ctxt2.getContext().getEA().decrypt(ctxt2, secret_key, ptxt_res2);
    for (int i = 0; i < nslots; i++)
        if (ptxt_res1[i] != ptxt_res2[i])
            return false;
    return true;
}

// Perform some operations before the actual bootstrapping to make it look real
bool squareWithThinBoot(helib::PubKey& pk, helib::Ctxt& c, const helib::SecKey* secret_key, bool our_version) {
    // Test polynomial evaluation with x^2 - x (for digit removal) and simultaneously
    // with f8 = 14641*x^8 + 22748*x^6 + 8836*x^4 + 112*x^2 (for digit extraction)
    /*****
    std::vector<NTL::ZZX> poly{NTL::ZZX(), NTL::ZZX()};
    SetCoeff(poly[0], 1, -1);
    SetCoeff(poly[0], 2, 1);
    SetCoeff(poly[1], 2, 112);
    SetCoeff(poly[1], 4, 8836);
    SetCoeff(poly[1], 6, 22748);
    SetCoeff(poly[1], 8, 14641);

    std::vector<helib::Ctxt> result;
    customPolyEval(result, poly, c);
    printPlaintext(result[0], *secret_key, "result 0");
    printPlaintext(result[1], *secret_key, "result 1");
    *****/

    bool refreshed = false;
    if (c.bitCapacity() <= 200) {
        helib::Ctxt tmp(c);
        pk.thinReCrypt(c, our_version, /*lazy*/ false, /*iterations*/ 5);
        if (areEqualCiphertexts(tmp, c, *secret_key))
            std::cout << "Thin bootstrapping successful!" << std::endl;
        else
            std::cout << "Thin bootstrapping failure!" << std::endl;
        refreshed = true;
    }

    helib::EncryptedArray ea(c.getContext().getEA()); int nslots = ea.size();
    std::vector<long> ptxt(nslots);
    for (int i = 0; i < nslots; ++i) {
        ptxt[i] = std::rand() % 256; // Random integers in plaintext space
    }
    helib::EncodedPtxt encoded_ptxt; ea.encode(encoded_ptxt, ptxt);
    c.square(); c.addConstant(encoded_ptxt);
    return refreshed;
}

// Perform some operations before the actual bootstrapping to make it look real
bool squareWithFatBoot(helib::PubKey& pk, helib::Ctxt& c, const helib::SecKey* secret_key, bool our_version) {
    bool refreshed = false;
    if (c.bitCapacity() <= 50) {
        helib::Ctxt tmp(c);
        pk.reCrypt(c, our_version, /*lazy*/ false);
        if (areEqualCiphertexts(tmp, c, *secret_key))
            std::cout << "Fat bootstrapping successful!" << std::endl;
        else
            std::cout << "Fat bootstrapping failure!" << std::endl;
        refreshed = true;
    }

    helib::EncryptedArray ea(c.getContext().getEA()); int nslots = ea.size();
    std::vector<long> ptxt(nslots);
    for (int i = 0; i < nslots; ++i) {
        ptxt[i] = std::rand() % 256; // Random integers in plaintext space
    }
    helib::EncodedPtxt encoded_ptxt; ea.encode(encoded_ptxt, ptxt);
    c.square(); c.addConstant(encoded_ptxt);
    return refreshed;
}

// Test bootstrapping with the given parameters
static void BM_thinboot(long m,
                        long p,
                        long r,
                        long c,
                        long bits,
                        long t,
                        int c_m,
                        std::vector<long> mvec,
                        std::vector<long> gens,
                        std::vector<long> ords,
                        bool our_version = true) {
    // Clang-format off
    std::cout << "m=" << m
            << ", p=" << p
            << ", e=" << r
            << std::endl;

    // Clang-format on
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .gens(gens)
                               .ords(ords)
                               .bits(bits)
                               .c(c)
                               .bootstrappable(true)
                               .skHwt(t)
                               .mvec(mvec)
                               .build();

    // Print the context
    if (our_version)
        std::cout << "Thin bootstrapping with our version" << std::endl;
    else
        std::cout << "Thin bootstrapping with built-in version" << std::endl;
    //std::cout << "Security: " << context.securityLevel() << std::endl;
    helib::SecKey secret_key(context);
    secret_key.GenSecKey();
    addSome1DMatrices(secret_key);
    addFrbMatrices(secret_key);
    // Generate bootstrapping data
    secret_key.genRecryptData();

    // NOTE: For some reason the reCrypt method is not marked const so
    //       I had to remove the const from the public key
    helib::PubKey& public_key = secret_key;
    const helib::EncryptedArray& ea = context.getEA();

    long nslots = ea.size();
    //std::cout << "Number of slots: " << nslots << std::endl;

    std::vector<long> ptxt(nslots);
    for (int i = 0; i < nslots; ++i) {
        ptxt[i] = std::rand() % 256; // Random integers in plaintext space
    }

    helib::Ctxt ctxt(public_key);
    ea.encrypt(ctxt, public_key, ptxt);
    long it = 0;
    while (!squareWithThinBoot(public_key, ctxt, &secret_key, our_version))
        it++;
    std::cout << std::endl;
}

// Test bootstrapping with the given parameters
static void BM_fatboot(long m,
                       long p,
                       long r,
                       long c,
                       long bits,
                       long t,
                       int c_m,
                       std::vector<long> mvec,
                       std::vector<long> gens,
                       std::vector<long> ords,
                       bool our_version = true) {
    // clang-format off
    std::cout << "m=" << m
            << ", p=" << p
            << ", e=" << r
            << std::endl;

    // clang-format on
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .gens(gens)
                               .ords(ords)
                               .bits(bits)
                               .c(c)
                               .bootstrappable(true)
                               .skHwt(t)
                               .mvec(mvec)
                               .thickboot()
                               .build();

    // Print the context
    if (our_version)
        std::cout << "Fat bootstrapping with our version" << std::endl;
    else
        std::cout << "Fat bootstrapping with built-in version" << std::endl;
    //std::cout << "Security: " << context.securityLevel() << std::endl;
    helib::SecKey secret_key(context);
    secret_key.GenSecKey();
    addSome1DMatrices(secret_key);
    addFrbMatrices(secret_key);
    secret_key.genRecryptData();

    // NOTE: For some reason the reCrypt method is not marked const so
    //       I had to remove the const from the public key
    helib::PubKey& public_key = secret_key;
    const helib::EncryptedArray& ea = context.getEA();

    long nslots = ea.size();
    //std::cout << "Number of slots: " << nslots << std::endl;

    std::vector<long> ptxt(nslots);
    for (int i = 0; i < nslots; ++i) {
        ptxt[i] = std::rand() % 256; // Random 0s and 1s
    }

    helib::Ctxt ctxt(public_key);
    ea.encrypt(ctxt, public_key, ptxt);
    long it = 0;
    while (!squareWithFatBoot(public_key, ctxt, &secret_key, our_version))
        it++;
    std::cout << std::endl;
}

// Test digit extraction with the given parameters
static void BM_thinextract(long m,
                           long p,
                           long r,
                           long c,
                           long bits,
                           long t,
                           int c_m,
                           int botHigh,
                           bool our_version = true,
                           std::vector<std::vector<long>> e_inner_compose_list = {{1}}) {
    // Clang-format off
    std::cout << "m=" << m
            << ", p=" << p
            << ", e=" << r
            << std::endl;

    // Clang-format on
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .bits(bits)
                               .c(c)
                               .bootstrappable(false)
                               .skHwt(t)
                               .build();

    // Print the context
    if (our_version) {
        std::cout << "Digit extraction with our version and e_inner_compose_list:" << std::endl;
        for (auto list : e_inner_compose_list) {
            for (auto element : list)
                std::cout << element << " ";

            std::cout << std::endl;
        }
    } else {
        std::cout << "Digit extraction with built-in version" << std::endl;
    }
    //std::cout << "Security: " << context.securityLevel() << std::endl;
    helib::SecKey secret_key(context);
    secret_key.GenSecKey();

    // NOTE: For some reason the reCrypt method is not marked const so
    //       I had to remove the const from the public key
    helib::PubKey& public_key = secret_key;
    const helib::EncryptedArray& ea = context.getEA();

    long nslots = ea.size();
    //std::cout << "Number of slots: " << nslots << std::endl;

    std::vector<long> ptxt(nslots);
    for (int i = 0; i < nslots; ++i) {
        ptxt[i] = std::rand() % 256; // Random integers in plaintext space
    }

    helib::Ctxt ctxt(public_key);
    ea.encrypt(ctxt, public_key, ptxt);

    // Test digit extraction
    wrapExtractDigitsThin(ctxt, botHigh, /*original ptxt exponent*/ r - botHigh, our_version, /*lazy*/ false, e_inner_compose_list, /*iterations*/ 5);

    std::vector<long> ptxt_res(nslots);
    ea.decrypt(ctxt, secret_key, ptxt_res);
    bool result = true;
    for (int i = 0; i < nslots; i++)
        if (((ptxt[i] + ((long)std::pow(p, botHigh)) / 2) / ((long)std::pow(p, botHigh))) != ptxt_res[i])
            result = false;
    if (result)
        std::cout << "Digit extraction successful!" << std::endl;
    else
        std::cout << "Digit extraction failure!" << std::endl;
    std::cout << std::endl;
}



// Useful parameters can be found at https://github.com/homenc/HElib/blob/master/tests/GTestThinBootstrapping.cpp
int main() {
// TOY PARAMETERS (not included in paper)

    BM_thinboot(  /*m = */      105,
                  /*p = */      2,
                  /*r = */      20,
                  /*c = */      3,
                  /*bits = */   1200,
                  /*t = */      120,
                  /*c_m = */    100,
                  /*mvec = */   std::vector<long>{3, 35},
                  /*gens =*/    std::vector<long>{71, 76},
                  /*ords =*/    std::vector<long>{2, 2});

    BM_fatboot(   /*m = */      105,
                  /*p = */      2,
                  /*r = */      20,
                  /*c = */      3,
                  /*bits = */   1200,
                  /*t = */      120,
                  /*c_m = */    100,
                  /*mvec = */   std::vector<long>{3, 35},
                  /*gens =*/    std::vector<long>{71, 76},
                  /*ords =*/    std::vector<long>{2, 2});

    BM_thinextract(  /*m = */      105,
                     /*p = */      2,
                     /*r = */      59,
                     /*c = */      3,
                     /*bits = */   1200,
                     /*t = */      120,
                     /*c_m = */    100,
                     /*botHigh*/   8);

// TESTS FOR FAT BOOTSTRAPPING

//    BM_fatboot(   /*m = */      42799,
//                  /*p = */      2,
//                  /*r = */      8,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    200,
//                  /*mvec = */   std::vector<long>{127, 337},
//                  /*gens =*/    std::vector<long>{25276, 40133},
//                  /*ords =*/    std::vector<long>{126, 16},
//                  /*version*/   false);
//
//    BM_fatboot(   /*m = */      42799,
//                  /*p = */      2,
//                  /*r = */      8,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    200,
//                  /*mvec = */   std::vector<long>{127, 337},
//                  /*gens =*/    std::vector<long>{25276, 40133},
//                  /*ords =*/    std::vector<long>{126, 16},
//                  /*version*/   true);
//
//    BM_fatboot(   /*m = */      45551,
//                  /*p = */      17,
//                  /*r = */      4,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    200,
//                  /*mvec = */   std::vector<long>{101, 451},
//                  /*gens =*/    std::vector<long>{19394, 7677},
//                  /*ords =*/    std::vector<long>{100, 10},
//                  /*version*/   false);
//
//    BM_fatboot(   /*m = */      45551,
//                  /*p = */      17,
//                  /*r = */      4,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    200,
//                  /*mvec = */   std::vector<long>{101, 451},
//                  /*gens =*/    std::vector<long>{19394, 7677},
//                  /*ords =*/    std::vector<long>{100, 10},
//                  /*version*/   true);
//
//    BM_fatboot(   /*m = */      32551,
//                  /*p = */      127,
//                  /*r = */      2,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    100,
//                  /*mvec = */   std::vector<long>{43, 757},
//                  /*gens =*/    std::vector<long>{7571, 28768},
//                  /*ords =*/    std::vector<long>{42, 54},
//                  /*version*/   false);
//
//    BM_fatboot(   /*m = */      32551,
//                  /*p = */      127,
//                  /*r = */      2,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    100,
//                  /*mvec = */   std::vector<long>{43, 757},
//                  /*gens =*/    std::vector<long>{7571, 28768},
//                  /*ords =*/    std::vector<long>{42, 54},
//                  /*version*/   true);

// TESTS FOR THIN BOOTSTRAPPING

//    BM_thinboot(  /*m = */      42799,
//                  /*p = */      2,
//                  /*r = */      8,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    200,
//                  /*mvec = */   std::vector<long>{127, 337},
//                  /*gens =*/    std::vector<long>{25276, 40133},
//                  /*ords =*/    std::vector<long>{126, 16},
//                  /*version*/   false);
//
//    BM_thinboot(  /*m = */      42799,
//                  /*p = */      2,
//                  /*r = */      8,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    200,
//                  /*mvec = */   std::vector<long>{127, 337},
//                  /*gens =*/    std::vector<long>{25276, 40133},
//                  /*ords =*/    std::vector<long>{126, 16},
//                  /*version*/   true);
//
//    BM_thinboot(  /*m = */      45551,
//                  /*p = */      17,
//                  /*r = */      4,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    200,
//                  /*mvec = */   std::vector<long>{101, 451},
//                  /*gens =*/    std::vector<long>{19394, 7677},
//                  /*ords =*/    std::vector<long>{100, 10},
//                  /*version*/   false);
//
//    BM_thinboot(  /*m = */      45551,
//                  /*p = */      17,
//                  /*r = */      4,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    200,
//                  /*mvec = */   std::vector<long>{101, 451},
//                  /*gens =*/    std::vector<long>{19394, 7677},
//                  /*ords =*/    std::vector<long>{100, 10},
//                  /*version*/   true);
//
//    BM_thinboot(  /*m = */      32551,
//                  /*p = */      127,
//                  /*r = */      2,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    100,
//                  /*mvec = */   std::vector<long>{43, 757},
//                  /*gens =*/    std::vector<long>{7571, 28768},
//                  /*ords =*/    std::vector<long>{42, 54},
//                  /*version*/   false);
//
//    BM_thinboot(  /*m = */      32551,
//                  /*p = */      127,
//                  /*r = */      2,
//                  /*c = */      3,
//                  /*bits = */   1200,
//                  /*t = */      120,
//                  /*c_m = */    100,
//                  /*mvec = */   std::vector<long>{43, 757},
//                  /*gens =*/    std::vector<long>{7571, 28768},
//                  /*ords =*/    std::vector<long>{42, 54},
//                  /*version*/   true);

// TESTS FOR DIGIT EXTRACTION

//    BM_thinextract(  /*m = */      42799,
//                     /*p = */      2,
//                     /*r = */      59,
//                     /*c = */      3,
//                     /*bits = */   1200,
//                     /*t = */      120,
//                     /*c_m = */    200,
//                     /*botHigh*/   8,
//                     /*version*/   false);
//
//    BM_thinextract(  /*m = */      42799,
//                     /*p = */      2,
//                     /*r = */      59,
//                     /*c = */      3,
//                     /*bits = */   1200,
//                     /*t = */      120,
//                     /*c_m = */    200,
//                     /*botHigh*/   8,
//                     /*version*/   true, {{1}});
//
//    BM_thinextract(  /*m = */      42799,
//                     /*p = */      2,
//                     /*r = */      59,
//                     /*c = */      3,
//                     /*bits = */   1200,
//                     /*t = */      120,
//                     /*c_m = */    200,
//                     /*botHigh*/   8,
//                     /*version*/   true, {{1, 16}});
//
//    BM_thinextract(  /*m = */      42799,
//                     /*p = */      2,
//                     /*r = */      59,
//                     /*c = */      3,
//                     /*bits = */   1200,
//                     /*t = */      120,
//                     /*c_m = */    200,
//                     /*botHigh*/   8,
//                     /*version*/   true, {{1, 16}, {1, 16}, {1, 16}, {1, 16}, {1, 16}, {1, 16}, {1, 16}, {1}});
//
//    BM_thinextract(  /*m = */      63973,
//                     /*p = */      3,
//                     /*r = */      37,
//                     /*c = */      3,
//                     /*bits = */   1400,
//                     /*t = */      120,
//                     /*c_m = */    200,
//                     /*botHigh*/   5,
//                     /*version*/   false);
//
//    BM_thinextract(  /*m = */      63973,
//                     /*p = */      3,
//                     /*r = */      37,
//                     /*c = */      3,
//                     /*bits = */   1400,
//                     /*t = */      120,
//                     /*c_m = */    200,
//                     /*botHigh*/   5,
//                     /*version*/   true, {{1}});
//
//    BM_thinextract(  /*m = */      63973,
//                     /*p = */      3,
//                     /*r = */      37,
//                     /*c = */      3,
//                     /*bits = */   1400,
//                     /*t = */      120,
//                     /*c_m = */    200,
//                     /*botHigh*/   5,
//                     /*version*/   true, {{1, 6}});
//
//    BM_thinextract(  /*m = */      63973,
//                     /*p = */      3,
//                     /*r = */      37,
//                     /*c = */      3,
//                     /*bits = */   1400,
//                     /*t = */      120,
//                     /*c_m = */    200,
//                     /*botHigh*/   5,
//                     /*version*/   true, {{1, 6}, {1, 6}, {1, 6}, {1, 6}, {1}});

    (void)BM_thinboot;
    (void)BM_fatboot;
    (void)BM_thinextract;

	return 0;
}
