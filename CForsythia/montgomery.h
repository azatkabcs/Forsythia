#pragma once
#include "parameters.h"
#include "gfp2.h"
#include <iostream>

class MontgomeryCurve {
private:
    void init();
protected:
    GFp2Element Ap24;
    GFp2Element C24;
    GFp2Element A;
    GFp2Element C;
    GFp2Element Am24;
    GFp2Element ap24;
public:
    MontgomeryCurve(const GFp2Element&, const GFp2Element&);
    MontgomeryCurve(mpz_t, const GFp2Element&);
    MontgomeryCurve(const GFp2Element&, mpz_t);
    MontgomeryCurve(mpz_t, mpz_t);
    friend std::ostream& operator<<(std::ostream &out, const MontgomeryCurve &elem);
    GFp2Element jinv() const;
};

class MontgomeryPoint : public MontgomeryCurve {
private:
    GFp2Element X;
    GFp2Element Z;
public:
    MontgomeryPoint(const GFp2Element& oX, const GFp2Element& oZ, const GFp2Element& oA, const GFp2Element& oC) : MontgomeryCurve(oA, oC), X(oX), Z(oZ) {}
    MontgomeryPoint(const GFp2Element& oX, const GFp2Element& oZ, const GFp2Element& oA, mpz_t oC) : MontgomeryCurve(oA, oC), X(oX), Z(oZ) {}
    MontgomeryPoint(const GFp2Element& oX, const GFp2Element& oZ, mpz_t oA, const GFp2Element& oC) : MontgomeryCurve(oA, oC), X(oX), Z(oZ) {}
    MontgomeryPoint(const GFp2Element& oX, const GFp2Element& oZ, mpz_t oA, mpz_t oC) : MontgomeryCurve(oA, oC), X(oX), Z(oZ) {}
    friend std::ostream& operator<<(std::ostream &out, const MontgomeryPoint &elem);
    MontgomeryPoint mul2() const;
    MontgomeryPoint mul3() const;
    MontgomeryPoint mul2e(mpz_t) const;
    MontgomeryPoint mul3e(mpz_t) const;
};