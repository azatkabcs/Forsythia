#include "montgomery.h"

void MontgomeryCurve::init() {
    mpz_t n2, n4;
    mpz_init_set_ui(n2, 2);
    mpz_init_set_ui(n4, 4);
    //Ap24 = A + C*n2;
    //Am24 = A - C*n2;
    Ap24 = C * n2;
    Am24 = A;
    Am24 = Am24 - Ap24;
    Ap24 = Ap24 + A;
    C24 = C*n4;
    ap24 = Ap24 / C;
    mpz_clears(n2, n4, NULL);
}

MontgomeryCurve::MontgomeryCurve(const GFp2Element &otherA, const GFp2Element &otherC): A(otherA), C(otherC) {
    init();
}

MontgomeryCurve::MontgomeryCurve(mpz_t otherA, const GFp2Element &otherC): A(otherA), C(otherC) {
    init();
}

MontgomeryCurve::MontgomeryCurve(const GFp2Element &otherA, mpz_t otherC): A(otherA), C(otherC) {
    init();
}

MontgomeryCurve::MontgomeryCurve(mpz_t otherA, mpz_t otherC): A(otherA), C(otherC) {
    init();
}

std::ostream& operator<<(std::ostream &out, const MontgomeryCurve &elem) {
    out << "y^2 = x^3 + (" << elem.A / elem.C << ") * x^2 + x";
    return out;
}

/**
 *  Get a Montgomery curve's j-invariant
 *  Alg. 9 from [SIKE]
 *  return: j-invariant of the curve (GFp2element)
 */
GFp2Element MontgomeryCurve::jinv() const{
    GFp2Element j = A * A;
    GFp2Element t1 = C * C;
    GFp2Element t0 = t1 + t1;
    t0 = j - t0;
    t0 = t0 - t1;
    j = t0 - t1;
    t1 = t1 * t1;
    j = j * t1;
    t0 = t0 + t0;
    t0 = t0 + t0;
    t1 = t0 * t0;
    t0 = t1 * t0;
    t0 = t0 + t0;
    t0 = t0 + t0;
    j = j.modinv();
    j = t0 * j;
    return GFp2Element(j);
}

std::ostream& operator<<(std::ostream& out, const MontgomeryPoint &elem) {
    out << "(" << elem.X << ":" << elem.Z << "); x = " << elem.X / elem.Z;
    return out;
}

/**
 *  Montgomery point x-only multiplication by 2
 *  Alg. 3 from [SIKE]
 */
MontgomeryPoint MontgomeryPoint::mul2() const {
    GFp2Element t0 = X - Z;
    GFp2Element t1 = X + Z;
    t0 = t0 * t0;
    t1 = t1 * t1;
    GFp2Element z = t0 * C24;
    GFp2Element x = z * t1;
    t1 = t1 - t0;
    t0 = t1 * Ap24;
    z = z + t0;
    z = z * t1;
    return MontgomeryPoint(x, z, A, C);
}

/**
 *  Montgomery point x-only multiplication by 3
 *  Alg. 6 from [SIKE]
 */
MontgomeryPoint MontgomeryPoint::mul3() const {
    GFp2Element t0 = X - Z;
    GFp2Element t2 = t0 * t0;
    GFp2Element t1 = X + Z;
    GFp2Element t3 = t1 * t1;
    GFp2Element t4 = t1 + t0;
    t0 = t1 - t0;
    t1 = t4 * t4;
    t1 = t1 - t3;
    t1 = t1 - t2;
    GFp2Element t5 = t3 * Ap24;
    t3 = t5 * t3;
    GFp2Element t6 = t2 * Am24;
    t2 = t2 * t6;
    t3 = t2 - t3;
    t2 = t5 - t6;
    t1 = t2 * t1;
    t2 = t3 + t1;
    t2 = t2 * t2;
    GFp2Element x = t2 * t4;
    t1 = t3 - t1;
    t1 = t1 * t1;
    GFp2Element z = t1 * t0;
    return MontgomeryPoint(x, z, A, C);
}

/**
 *  e-repeated Montgomery point x-only multiplication by 2
 *  Alg. 4 from [SIKE]
 */
MontgomeryPoint MontgomeryPoint::mul2e(mpz_t e) const {
    MontgomeryPoint res(X, Z, A, C);
    mpz_t index;
    for (mpz_init(index); mpz_cmp(index, e) < 0; mpz_add_ui(index, index, 1))
        res = res.mul2();
    mpz_clear(index);
    return res;
}

/**
 *  e-repeated Montgomery point x-only multiplication by 3
 *  Alg. 7 from [SIKE]
 */
MontgomeryPoint MontgomeryPoint::mul3e(mpz_t e) const {
    MontgomeryPoint res(X, Z, A, C);
    mpz_t index;
    for (mpz_init(index); mpz_cmp(index, e) < 0; mpz_add_ui(index, index, 1))
        res = res.mul3();
    mpz_clear(index);
    return res;
}