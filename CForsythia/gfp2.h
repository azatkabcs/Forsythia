#pragma once
#include <gmp.h>
#include "parameters.h"
#include <iostream>

class GFp2Element {
private:
    mpz_t a;
    mpz_t b;
public:
    GFp2Element();
    GFp2Element(mpz_t a, mpz_t b);
    GFp2Element(mpz_t a);
    GFp2Element(unsigned int);
    GFp2Element(const GFp2Element&);
    ~GFp2Element();
    GFp2Element operator+(const GFp2Element&) const;
    GFp2Element operator+(mpz_t);
    GFp2Element operator-(const GFp2Element&) const;
    GFp2Element operator-(mpz_t);
    GFp2Element operator*(const GFp2Element&) const;
    GFp2Element operator*(mpz_t);
    GFp2Element modinv() const;
    GFp2Element operator/(const GFp2Element&) const;
    GFp2Element operator/(mpz_t);
    GFp2Element& operator=(const GFp2Element&);
    bool operator==(const GFp2Element&) const;
    bool operator==(mpz_t);
    bool operator!=(const GFp2Element&) const;
    bool operator!=(mpz_t);
    bool iszero() const;
    friend std::ostream& operator<<(std::ostream& out, const GFp2Element& elem);
};