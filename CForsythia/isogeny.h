#pragma once
#include "montgomery.h"

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
);

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
);

GFp2Element isoex2(mpz_t sk2, mpz_t e2, GFp2Element pk[3]);
GFp2Element isoex3(mpz_t sk3, mpz_t e3, GFp2Element pk[3]);