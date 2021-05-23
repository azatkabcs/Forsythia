#pragma once
#include "parameters.h"
#include "gfp2.h"
#include <iostream>

class MontgomeryPoint;

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
    MontgomeryCurve() {}
    MontgomeryCurve(const GFp2Element&, const GFp2Element&);
    MontgomeryCurve(mpz_t, const GFp2Element&);
    MontgomeryCurve(const GFp2Element&, mpz_t);
    MontgomeryCurve(mpz_t, mpz_t);
    MontgomeryCurve(const MontgomeryCurve&);
    friend std::ostream& operator<<(std::ostream &out, const MontgomeryCurve &elem);
    GFp2Element jinv() const;
    void xdbladd(const MontgomeryPoint&, const MontgomeryPoint&, const MontgomeryPoint&, MontgomeryPoint&, MontgomeryPoint&) const;
    void ladder3pt(mpz_t, const GFp2Element&, const GFp2Element&, const GFp2Element&, MontgomeryPoint&);
    void seta(const GFp2Element&, const GFp2Element&, const GFp2Element&);
    MontgomeryCurve iso2_curve(const MontgomeryPoint&) const;
    MontgomeryPoint iso2_eval(const MontgomeryPoint&, const MontgomeryPoint&, const MontgomeryCurve&) const;
    void iso4_curve(const MontgomeryPoint&, MontgomeryCurve&, GFp2Element&, GFp2Element&, GFp2Element&) const;
    MontgomeryPoint iso4_eval(const GFp2Element&, const GFp2Element&, const GFp2Element&, const MontgomeryPoint&, const MontgomeryCurve&) const;
    void iso3_curve(const MontgomeryPoint&, MontgomeryCurve&, GFp2Element&, GFp2Element&) const;
    MontgomeryPoint iso3_eval(const GFp2Element&, const GFp2Element&, const MontgomeryPoint&, const MontgomeryCurve&) const;
    void iso2e(mpz_t, const MontgomeryPoint&, const GFp2Element&, const GFp2Element&, const GFp2Element&, MontgomeryCurve&, MontgomeryPoint&, MontgomeryPoint&, MontgomeryPoint&) const;
    void iso2eby4(mpz_t, const MontgomeryPoint&, const GFp2Element&, const GFp2Element&, const GFp2Element&, MontgomeryCurve&, MontgomeryPoint&, MontgomeryPoint&, MontgomeryPoint&) const;
    void iso3e(mpz_t, const MontgomeryPoint&, const GFp2Element&, const GFp2Element&, const GFp2Element&, MontgomeryCurve&, MontgomeryPoint&, MontgomeryPoint&, MontgomeryPoint&) const;
};

class MontgomeryPoint : public MontgomeryCurve {
private:
    GFp2Element X;
    GFp2Element Z;
public:
    MontgomeryPoint() {}
    MontgomeryPoint(const MontgomeryPoint& other): MontgomeryCurve(other.A, other.C), X(other.X), Z(other.Z) {}
    MontgomeryPoint(const GFp2Element& oX, const GFp2Element& oZ, const GFp2Element& oA, const GFp2Element& oC) : MontgomeryCurve(oA, oC), X(oX), Z(oZ) {}
    MontgomeryPoint(const GFp2Element& oX, const GFp2Element& oZ, const GFp2Element& oA, mpz_t oC) : MontgomeryCurve(oA, oC), X(oX), Z(oZ) {}
    MontgomeryPoint(const GFp2Element& oX, const GFp2Element& oZ, mpz_t oA, const GFp2Element& oC) : MontgomeryCurve(oA, oC), X(oX), Z(oZ) {}
    MontgomeryPoint(const GFp2Element& oX, const GFp2Element& oZ, mpz_t oA, mpz_t oC) : MontgomeryCurve(oA, oC), X(oX), Z(oZ) {}
    MontgomeryPoint(const GFp2Element& oX, const GFp2Element& oZ, const MontgomeryCurve& curve) : MontgomeryCurve(curve), X(oX), Z(oZ) {}
    
    friend std::ostream& operator<<(std::ostream &out, const MontgomeryPoint &elem);
    MontgomeryPoint mul2() const;
    MontgomeryPoint mul3() const;
    MontgomeryPoint mul2e(mpz_t) const;
    MontgomeryPoint mul3e(mpz_t) const;
    const GFp2Element getX() const { return X; }
    const GFp2Element getZ() const { return Z; }
    GFp2Element getNatureX() const { return X / Z; }
};
