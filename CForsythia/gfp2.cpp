#include "gfp2.h"

GFp2Element::GFp2Element(mpz_t x, mpz_t y) {
    mpz_inits(a, b, NULL);
    mpz_set(a, x);
    mpz_set(b, y);
}

GFp2Element GFp2Element::operator+(GFp2Element &other) {
    mpz_t buffer_a, buffer_b;
    mpz_inits(buffer_a, buffer_b, NULL);
    mpz_add(buffer_a, this->a, other.a);
    mpz_mod(buffer_a, buffer_a, p_global);
    mpz_add(buffer_b, this->b, other.b);
    mpz_mod(buffer_b, buffer_b, p_global);
    return GFp2Element(buffer_a, buffer_b);
}

GFp2Element GFp2Element::operator+(mpz_t other) {
    mpz_t buffer_a;
    mpz_init(buffer_a);
    mpz_add(buffer_a, this->a, other);
    mpz_mod(buffer_a, buffer_a, p_global);
    return GFp2Element(buffer_a, this->b);
}

GFp2Element GFp2Element::operator-(GFp2Element &other) {
    mpz_t buffer_a, buffer_b;
    mpz_inits(buffer_a, buffer_b, NULL);
    mpz_sub(buffer_a, this->a, other.a);
    mpz_mod(buffer_a, buffer_a, p_global);
    mpz_sub(buffer_b, this->b, other.b);
    mpz_mod(buffer_b, buffer_b, p_global);
    return GFp2Element(buffer_a, buffer_b);
}

GFp2Element GFp2Element::operator-(mpz_t other) {
    mpz_t buffer_a;
    mpz_init(buffer_a);
    mpz_sub(buffer_a, this->a, other);
    mpz_mod(buffer_a, buffer_a, p_global);
    return GFp2Element(buffer_a, this->b);
}

GFp2Element GFp2Element::operator*(GFp2Element &other) {
    mpz_t buffer_a, buffer_b, tmp1, tmp2;
    mpz_inits(buffer_a, buffer_b, tmp1, tmp2, NULL);
    mpz_mul(buffer_a, this->a, other.a);
    mpz_mul(buffer_b, this->b, other.b);
    mpz_sub(tmp1, buffer_a, buffer_b);
    mpz_mul(buffer_a, this->a, other.b);
    mpz_mul(buffer_b, this->b, other.a);
    mpz_add(tmp2, buffer_a, buffer_b);
    mpz_mod(buffer_a, tmp1, p_global);
    mpz_mod(buffer_b, tmp2, p_global);
    return GFp2Element(buffer_a, buffer_b);
}

GFp2Element GFp2Element::operator*(mpz_t other) {
    mpz_t buffer_a, buffer_b;
    mpz_inits(buffer_a, buffer_b, NULL);
    mpz_mul(buffer_a, this->a, other);
    mpz_mod(buffer_a, buffer_a, p_global);
    mpz_mul(buffer_b, this->b, other);
    mpz_mod(buffer_b, buffer_b, p_global);
    return GFp2Element(buffer_a, buffer_b);
}

GFp2Element GFp2Element::modinv() {
    mpz_t buffer_a, buffer_b, buffer, j;
    mpz_inits(buffer_a, buffer_b, buffer, j, NULL);
    mpz_mul(buffer, this->a, this->a);
    mpz_mul(j, this->b, this->b);
    mpz_add(buffer, buffer, j);
    mpz_mod(buffer, buffer, p_global);
    mpz_invert(j, buffer, p_global);
    mpz_mul(buffer_a, this->a, j);
    mpz_mod(buffer_a, buffer_a, p_global);
    mpz_neg(buffer_b, this->b);
    mpz_mul(buffer_b, buffer_b, j);
    mpz_mod(buffer_b, buffer_b, p_global);
    return GFp2Element(buffer_a, buffer_b);
}

GFp2Element GFp2Element::operator/(GFp2Element &other) {
    GFp2Element buffer = other.modinv();
    return this->operator*(buffer);
}

GFp2Element GFp2Element::operator/(mpz_t other) {
    mpz_t buffer;
    mpz_init_set(buffer, other);
    mpz_invert(buffer, buffer, p_global);
    return this->operator*(buffer);
}

bool GFp2Element::operator==(GFp2Element &other) {
    return !mpz_cmp(a, other.a) && !mpz_cmp(b, other.b);
}

bool GFp2Element::operator==(mpz_t other) {
    return !mpz_cmp(a, other) && !mpz_cmp_ui(b, 0);
}

bool GFp2Element::operator!=(GFp2Element &other) {
    return !this->operator==(other);
}

bool GFp2Element::operator!=(mpz_t other) {
    return !this->operator==(other);
}

bool GFp2Element::iszero() {
    return !mpz_cmp_ui(a, 0) && !mpz_cmp_ui(b, 0);
}

std::ostream& operator<<(std::ostream& out, const GFp2Element& elem) {
    gmp_printf("a: %Zd\n", elem.a);
    gmp_printf("b: %Zd", elem.b);
    return out;
}