#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <gmp.h>

using namespace std;

void buildsmallprimes(uint32_t c, vector<uint32_t>& result);
bool MR(mpz_t c, uint32_t k);