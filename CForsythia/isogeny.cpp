#include "isogeny.h"

/**
 *  Generate public key in 2^e-torsion
 *  Alg. 21 from [SIKE]
 *  :param e0: starting curve
 *  :param sk2: Alice's secret key
 *  :param e2: Degree of 2
 *  :param xp2: X-coordinate of Alice's basis point P
 *  :param xq2: X-coordinate of Alice's basis point Q
 *  :param xr2: X-coordinate of Alice's basis point Q-P
 *  :param xp3: X-coordinate of Bob's basis point P
 *  :param xq3: X-coordinate of Bob's basis point Q
 *  :param xr3: X-coordinate of Bob's basis point Q-P
 */
void isogen2(
    MontgomeryCurve& e0,
    mpz_t sk2,
    mpz_t e2,
    const GFp2Element& xp2,
    const GFp2Element& xq2,
    const GFp2Element& xr2,
    const GFp2Element& xp3,
    const GFp2Element& xq3,
    const GFp2Element& xr3,
    GFp2Element& x1,
    GFp2Element& x2,
    GFp2Element& x3
) {
    MontgomeryPoint S, X1, X2, X3;
    MontgomeryCurve curve;
    e0.ladder3pt(sk2, xp2, xq2, xr2, S);
    e0.iso2e(e2, S, xp3, xq3, xr3, curve, X1, X2, X3);
    x1 = X1.getNatureX();
    x2 = X2.getNatureX();
    x3 = X3.getNatureX();
}

/**
 *  Generate public key in 3^e3-torsion
 *  Alg. 22 from [SIKE]
 *  :param e0: starting curve
 *  :param sk3: Bob's secret key
 *  :param e3: Degree of 3
 *  :param xp2: X-coordinate of Alice's basis point P
 *  :param xq2: X-coordinate of Alice's basis point Q
 *  :param xr2: X-coordinate of Alice's basis point Q-P
 *  :param xp3: X-coordinate of Bob's basis point P
 *  :param xq3: X-coordinate of Bob's basis point Q
 *  :param xr3: X-coordinate of Bob's basis point Q-P
 */
void isogen3(
    MontgomeryCurve& e0,
    mpz_t sk3,
    mpz_t e3,
    const GFp2Element& xp2,
    const GFp2Element& xq2,
    const GFp2Element& xr2,
    const GFp2Element& xp3,
    const GFp2Element& xq3,
    const GFp2Element& xr3,
    GFp2Element& x1,
    GFp2Element& x2,
    GFp2Element& x3
) {
    MontgomeryPoint S, X1, X2, X3;
    MontgomeryCurve curve;
    e0.ladder3pt(sk3, xp3, xq3, xr3, S);
    e0.iso3e(e3, S, xp2, xq2, xr2, curve, X1, X2, X3);
    x1 = X1.getNatureX();
    x2 = X2.getNatureX();
    x3 = X3.getNatureX();
}

/**
 *  Generate shared key in 2^e2-torsion
 *  Alg. 23 from [SIKE]
 *  :param sk2: Alice's secret key
 *  :param e2: Power of 2
 *  :param pk: Bob's public key encoded as three points
 */
GFp2Element isoex2(mpz_t sk2, mpz_t e2, GFp2Element pk[3]) {
    mpz_t unit;
    mpz_init_set_ui(unit, 1);
    MontgomeryCurve curve(unit, unit);
    mpz_clear(unit);
    GFp2Element x1 = pk[0], x2 = pk[1], x3 = pk[2];
    curve.seta(x1, x2, x3);
    MontgomeryPoint S;
    curve.ladder3pt(sk2, x1, x2, x3, S);
    MontgomeryCurve image;
    GFp2Element buf1, buf2, buf3;
    MontgomeryPoint X1, X2, X3;
    curve.iso2e(e2, S, buf1, buf2, buf3, image, X1, X2, X3);

    return image.jinv();
}

/**
 *  Generate shared key in 3^e3-torsion
 *  Alg. 24 from [SIKE]
 *  :param sk3: Bob's secret key
 *  :param e3: Power of 3
 *  :param pk: Alice's public key encoded as three points
 */
GFp2Element isoex3(mpz_t sk3, mpz_t e3, GFp2Element pk[3]) {
    mpz_t unit;
    mpz_init_set_ui(unit, 1);
    MontgomeryCurve curve(unit, unit);
    mpz_clear(unit);
    GFp2Element x1 = pk[0], x2 = pk[1], x3 = pk[2];
    curve.seta(x1, x2, x3);
    MontgomeryPoint S;
    curve.ladder3pt(sk3, x1, x2, x3, S);
    MontgomeryCurve image;
    GFp2Element buf1, buf2, buf3;
    MontgomeryPoint X1, X2, X3;
    curve.iso3e(e3, S, buf1, buf2, buf3, image, X1, X2, X3);

    return image.jinv();
}