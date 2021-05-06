#pragma once
#include <gmp.h>
#include "parameters.h"

class GFp2Element {
private:
    mpz_t a;
    mpz_t b;
public:
    GFp2Element(mpz_t a, mpz_t b);
    GFp2Element operator+(GFp2Element&);
    GFp2Element operator+(mpz_t);
    GFp2Element operator-(GFp2Element&);
    GFp2Element operator-(mpz_t);
    GFp2Element operator*(GFp2Element&);
    GFp2Element operator*(mpz_t);
    GFp2Element modinv();
    GFp2Element operator/(GFp2Element&);
    GFp2Element operator/(mpz_t);
    bool operator==(GFp2Element&);
    bool operator==(mpz_t);
    bool operator!=(GFp2Element&);
    bool operator!=(mpz_t);
    bool iszero();
    friend std::ostream& operator<<(std::ostream& out, const GFp2Element& elem);
};