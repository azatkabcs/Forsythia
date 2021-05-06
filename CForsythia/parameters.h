#pragma once
#include <stdint.h>
#include <gmp.h>

typedef struct point_t {
    const char *first;
    const char *second;
} point;

typedef struct parameters_t {
    uint32_t eA;
    uint32_t eB;
    uint32_t f;
    point A;
    point C;
    point xp2;
    point yp2;
    point zp2;
    point xq2;
    point yq2;
    point zq2;
    point xp3;
    point yp3;
    point zp3;
    point xq3;
    point yq3;
    point zq3;
} parameters;

extern mpz_t p_global;
extern parameters forsythia80;
