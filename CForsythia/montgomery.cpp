#include "montgomery.h"

void MontgomeryCurve::init() {
    mpz_t n2, n4;
    mpz_init_set_ui(n2, 2);
    mpz_init_set_ui(n4, 4);
    Ap24 = A + C*n2;
    Am24 = A - C*n2;
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

MontgomeryCurve::MontgomeryCurve(const MontgomeryCurve& other): A(other.A), C(other.C), Ap24(other.Ap24), C24(other.C24), Am24(other.Am24), ap24(other.ap24) {}

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

/**
 * Double-and-Add, Alg. 5 in [SIKE]
 */
void MontgomeryCurve::xdbladd(const MontgomeryPoint &P, const MontgomeryPoint &Q, const MontgomeryPoint &D, MontgomeryPoint& p0, MontgomeryPoint& p1) const {
    GFp2Element t0(P.getX() + P.getZ());
    GFp2Element t1(P.getX() - P.getZ());
    GFp2Element x2p(t0 * t0);
    GFp2Element t2(Q.getX() - Q.getZ());
    GFp2Element xpq(Q.getX() + Q.getZ());
    t0 = t0 * t2;
    GFp2Element z2p(t1 * t1);
    t1 = t1 * xpq;
    t2 = x2p - z2p;
    x2p = x2p * z2p;
    xpq = ap24 * t2;
    GFp2Element zpq(t0 - t1);
    z2p = xpq + z2p;
    xpq = t0 + t1;
    z2p = z2p * t2;
    zpq = zpq * zpq;
    xpq = xpq * xpq;
    zpq = D.getX() * zpq;
    xpq = D.getZ() * xpq;
    p0 = MontgomeryPoint(x2p, z2p, A, C);
    p1 = MontgomeryPoint(xpq, zpq, A, C);
}

/**
 *  Montgomery's ladder. Calculates x(P+[m]Q) given m and x-coordinates of P, Q, D=Q-P
 *  Alg. 9 in [SIKE]
 */
void MontgomeryCurve::ladder3pt(mpz_t m, const GFp2Element& xP, const GFp2Element& xQ, const GFp2Element& xD, MontgomeryPoint &res) {
    mpz_t unit, n2, n4, index;
    mpz_init_set_ui(unit, 1);
    mpz_init_set_ui(n2, 2);
    mpz_init_set_ui(n4, 4);
    mpz_init_set(index, m);
    GFp2Element elem_unit(unit);
    MontgomeryPoint p0(xQ, elem_unit, A, C);
    MontgomeryPoint p1(xP, elem_unit, A, C);
    MontgomeryPoint p2(xD, elem_unit, A, C);
    MontgomeryPoint res0, res1;
    ap24 = A + n2;
    ap24 = ap24 / n4;
    while (mpz_sgn(index))
    {
        if (mpz_fdiv_ui(index, 2)) {
            xdbladd(p0, p1, p2, res0, res1);
            p0 = res0;
            p1 = res1;
        }
        else {
            xdbladd(p0, p2, p1, res0, res1);
            p0 = res0;
            p2 = res1;
        }
        mpz_fdiv_q_2exp(index, index, 1);
    }
    mpz_clears(unit, n2, n4, index, NULL);
    res = p1;
}

/**
 *  Recover Montgomery curve coefficient A as well as aux curve constants from P, Q, P-Q x-coordinates
 *  Alg. 10 in [SIKE]
 */
void MontgomeryCurve::seta(const GFp2Element& p, const GFp2Element& q, const GFp2Element& d) {
    GFp2Element t1(p + q);
    GFp2Element t0(p * q);
    GFp2Element _A(d * t1);
    mpz_t unit, n2, n4;
    mpz_init_set_ui(unit, 1);
    mpz_init_set_ui(n2, 2);
    mpz_init_set_ui(n4, 4);
    _A = _A + t0;
    t0 = t0 * d;
    _A = _A - unit;
    t0 = t0 + t0;
    t1 = t1 + d;
    t0 = t0 + t0;
    _A = _A * _A;
    t0 = t0.modinv();
    _A = _A * t0;
    _A = _A - t1;
    A = _A;
    C = GFp2Element(unit);
    Ap24 = A + C * n2;
    Am24 = A - C * n2;
    C24 = C * n4;
    ap24 = Ap24 / C;
    mpz_clears(unit, n2, n4, NULL);
}

/**
 *  Calculate 2-isogenous curve
 *  Alg. 11 from [SIKE]
 */
MontgomeryCurve MontgomeryCurve::iso2_curve(const MontgomeryPoint &P2) const {
    mpz_t n2, n4;
    mpz_init_set_ui(n2, 2);
    mpz_init_set_ui(n4, 4);
    GFp2Element _Ap24(P2.getX() * P2.getX());
    GFp2Element _C24(P2.getZ() * P2.getZ());
    _Ap24 = _C24 - _Ap24;
    GFp2Element _A(_Ap24 * n4 - _C24 * n2);
    mpz_clears(n2, n4, NULL);
    return MontgomeryCurve(_A, _C24);
}

/**
 *  Evaluate a 2-isogeny on a point
 *  Alg. 12 from [SIKE]
 */
MontgomeryPoint MontgomeryCurve::iso2_eval(const MontgomeryPoint& P2, const MontgomeryPoint& Q, const MontgomeryCurve& image) const {
    GFp2Element t0(P2.getX() + P2.getZ());
    GFp2Element t1(P2.getX() - P2.getZ());
    GFp2Element t2(Q.getX() + Q.getZ());
    GFp2Element t3(Q.getX() - Q.getZ());
    t0 = t0 * t3;
    t1 = t1 * t2;
    t2 = t0 + t1;
    t3 = t0 - t1;
    GFp2Element XQP(Q.getX() * t2);
    GFp2Element ZQP(Q.getZ() * t3);
    return MontgomeryPoint(XQP, ZQP, image);
}

/**
 *  Calculate 4-isogenous curve
 *  Alg. 13 from [SIKE]
 */
void MontgomeryCurve::iso4_curve(const MontgomeryPoint &P4, MontgomeryCurve &curve, GFp2Element &K1, GFp2Element &K2, GFp2Element &K3) const {
    K2 = P4.getX() - P4.getZ();
    K3 = P4.getX() + P4.getZ();
    K1 = P4.getZ() * P4.getZ();
    K1 = K1 + K1;
    GFp2Element _C24(K1 * K1);
    K1 = K1 + K1;
    GFp2Element _Ap24(P4.getX() * P4.getX());
    _Ap24 = _Ap24 + _Ap24;
    _Ap24 = _Ap24 * _Ap24;
    mpz_t n2, n4;
    mpz_init_set_ui(n2, 2);
    mpz_init_set_ui(n4, 4);
    GFp2Element _A(_Ap24 * n4 - _C24 * n2);
    mpz_clears(n2, n4, NULL);
    curve = MontgomeryCurve(_A, _C24);
}

/**
 *  Evaluate a 4-isogeny at a point
 *  Alg. 14 from [SIKE] has a bug, don't know how to fix it =(
 */
MontgomeryPoint MontgomeryCurve::iso4_eval(const GFp2Element& K1, const GFp2Element& K2, const GFp2Element& K3, const MontgomeryPoint& Q, const MontgomeryCurve& image) const {
    GFp2Element t0(Q.getX() + Q.getZ());
    GFp2Element t1(Q.getX() - Q.getZ());
    GFp2Element XPQ(t0 * K2);
    GFp2Element ZPQ(t1 * K3);
    t0 = t0 * t1;
    t0 = t0 * K1;
    t1 = XPQ + ZPQ;
    ZPQ = XPQ - ZPQ;
    t1 = t1 * t1;
    ZPQ = ZPQ * ZPQ;
    XPQ = t0 + t1;
    t0 = ZPQ - t1;
    XPQ = XPQ * t1;
    ZPQ = ZPQ * t0;
    return MontgomeryPoint(XPQ, ZPQ, image);
}

/**
 *  Calculate 3-isogenous curve and parameters K1, K2
 *  Alg. 15 from [SIKE]
 */
void MontgomeryCurve::iso3_curve(const MontgomeryPoint& P3, MontgomeryCurve& curve, GFp2Element& K1, GFp2Element& K2) const {
    K1 = P3.getX() - P3.getZ();
    GFp2Element t0(K1 * K1);
    K2 = P3.getX() + P3.getZ();
    GFp2Element t1(K2 * K2);
    GFp2Element t2(t0 + t1);
    GFp2Element t3(K1 + K2);
    t3 = t3 * t3;
    t3 = t3 - t2;
    t2 = t1 + t3;
    t3 = t3 + t0;
    GFp2Element t4(t3 + t0);
    t4 = t4 + t4;
    t4 = t1 + t4;
    GFp2Element _Am24(t2 * t4);
    t4 = t1 + t2;
    t4 = t4 + t4;
    t4 = t0 + t4;
    GFp2Element _Ap24(t3 * t4);
    mpz_t n2;
    mpz_init_set_ui(n2, 2);
    GFp2Element _A(_Ap24 * n2 + _Am24 * n2);
    mpz_clear(n2);
    GFp2Element _C(_Ap24 - _Am24);
    curve = MontgomeryCurve(_A, _C);
}

/**
 *  Alg. 16 from [SIKE]
 */
MontgomeryPoint MontgomeryCurve::iso3_eval(const GFp2Element& K1, const GFp2Element& K2, const MontgomeryPoint& Q, const MontgomeryCurve& image) const {
    GFp2Element t0(Q.getX() + Q.getZ());
    GFp2Element t1(Q.getX() - Q.getZ());
    t0 = K1 * t0;
    t1 = K2 * t1;
    GFp2Element t2(t0 + t1);
    t0 = t1 - t0;
    t2 = t2 * t2;
    t0 = t0 * t0;
    GFp2Element XPQ(Q.getX() * t2);
    GFp2Element ZPQ(Q.getZ() * t0);
    return MontgomeryPoint(XPQ, ZPQ, image);
}

/**
 *  Compute and optionally evaluate a 2^e2-isogeny
 */
void MontgomeryCurve::iso2e(
    mpz_t e2, 
    const MontgomeryPoint& S1, 
    const GFp2Element& X11, 
    const GFp2Element& X22, 
    const GFp2Element& X33, 
    MontgomeryCurve& curve, 
    MontgomeryPoint& X1, 
    MontgomeryPoint& X2, 
    MontgomeryPoint& X3) const {
    
    MontgomeryPoint S(S1);
    mpz_t unit, index;
    mpz_init_set_ui(unit, 1);
    if (!X11.iszero()) {
        X1 = MontgomeryPoint(X11, GFp2Element(unit), *this);
    }
    if (!X22.iszero()) {
        X2 = MontgomeryPoint(X22, GFp2Element(unit), *this);
    }
    if (!X33.iszero()) {
        X3 = MontgomeryPoint(X33, GFp2Element(unit), *this);
    }
    MontgomeryPoint T;
    for (mpz_init_set(index, e2), mpz_sub_ui(index, index, 1); mpz_cmp_si(index, -1) > 0; mpz_sub_ui(index, index, 1)) {
        //gmp_printf("index: %Zd\n", index);
        T = S.mul2e(index);
        //std::cout << "T: " << T << std::endl;
        curve = iso2_curve(T);
        //std::cout << "curve: " << curve << std::endl;
        if (mpz_sgn(index)) {
            S = iso2_eval(T, S, curve);
            //std::cout << "S: "<< S << std::endl;
        }
        if (!X11.iszero()) {
            X1 = iso2_eval(T, X1, curve);
        }
        if (!X22.iszero()) {
            X2 = iso2_eval(T, X2, curve);
        }
        if (!X33.iszero()) {
            X3 = iso2_eval(T, X3, curve);
        }
    }
    mpz_clears(unit, index, NULL);
}

/**
 *  Compute and optionally evaluate a 2^e2-isogeny
 *  Alg. 17 from [SIKE]
 */
void MontgomeryCurve::iso2eby4(
    mpz_t e2, 
    const MontgomeryPoint& S, 
    const GFp2Element& X11, 
    const GFp2Element& X22, 
    const GFp2Element& X33, 
    MontgomeryCurve& curve, 
    MontgomeryPoint& X1, 
    MontgomeryPoint& X2, 
    MontgomeryPoint& X3) const {

    mpz_t unit, index;
    mpz_init_set_ui(unit, 1);
    if (!X11.iszero()) {
        X1 = MontgomeryPoint(X11, GFp2Element(unit), *this);
    }
    if (!X22.iszero()) {
        X2 = MontgomeryPoint(X22, GFp2Element(unit), *this);
    }
    if (!X33.iszero()) {
        X3 = MontgomeryPoint(X33, GFp2Element(unit), *this);
    }
    MontgomeryPoint T;
    MontgomeryPoint S1(S);
    GFp2Element K1, K2, K3;
    for (mpz_init_set(index, e2), mpz_sub_ui(index, index, 2); mpz_cmp_si(index, -2) > 0; mpz_sub_ui(index, index, 2)) {
        T = S1.mul2e(index);
        iso4_curve(T, curve, K1, K2, K3);
        if (mpz_sgn(index)) {
            S1 = iso4_eval(K1, K2, K3, S1, curve);
        }
        if (!X11.iszero()) {
            X1 = iso4_eval(K1, K2, K3, X1, curve);
        }
        if (!X22.iszero()) {
            X2 = iso4_eval(K1, K2, K3, X2, curve);
        }
        if (!X33.iszero()) {
            X3 = iso4_eval(K1, K2, K3, X3, curve);
        }
    }
    mpz_clears(unit, index, NULL);
}

/**
 *  Compute and optionally evaluate a 3^e-isogeny
 *  Alg. 18 from [SIKE]
 */
void MontgomeryCurve::iso3e(
    mpz_t e3, 
    const MontgomeryPoint& S1, 
    const GFp2Element& X11, 
    const GFp2Element& X22, 
    const GFp2Element& X33, 
    MontgomeryCurve& curve, 
    MontgomeryPoint& X1, 
    MontgomeryPoint& X2, 
    MontgomeryPoint& X3) const {

    MontgomeryPoint S(S1);
    mpz_t unit, index;
    mpz_init_set_ui(unit, 1);
    if (!X11.iszero()) {
        X1 = MontgomeryPoint(X11, GFp2Element(unit), *this);
    }
    if (!X22.iszero()) {
        X2 = MontgomeryPoint(X22, GFp2Element(unit), *this);
    }
    if (!X33.iszero()) {
        X3 = MontgomeryPoint(X33, GFp2Element(unit), *this);
    }
    MontgomeryPoint T;
    GFp2Element K1, K2;
    for (mpz_init_set(index, e3), mpz_sub_ui(index, index, 1); mpz_cmp_si(index, -1) > 0; mpz_sub_ui(index, index, 1)) {
        T = S.mul3e(index);
        iso3_curve(T, curve, K1, K2);
        if (mpz_sgn(index)) {
            S = iso3_eval(K1, K2, S, curve);
        }
        if (!X11.iszero()) {
            X1 = iso3_eval(K1, K2, X1, curve);
        }
        if (!X22.iszero()) {
            X2 = iso3_eval(K1, K2, X2, curve);
        }
        if (!X33.iszero()) {
            X3 = iso3_eval(K1, K2, X3, curve);
        }
    }
    mpz_clears(unit, index, NULL);
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