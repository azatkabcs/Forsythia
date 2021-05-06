#include "utils.h"

void buildsmallprimes(uint32_t c, vector<uint32_t> &result) {
	vector<bool> sprimes(c);

	for (uint32_t i = 2; i < c; i++)
		if (sprimes[i] == false)
			for (int j = 2 * i; j < c; j += i)
				sprimes[j] = true;

	for (uint32_t i = 2; i < c; i++)
		if (sprimes[i] == false)
			result.push_back(i);
}

/* Miller-Rabin */
bool MR(mpz_t c, uint32_t k) {
	mpz_t N, u, t, a, gcd, i, a_mod_N, N_minus_1;
	mpz_inits(N, u, t, a, gcd, i, a_mod_N, N_minus_1, NULL);
	mpz_set(N, c);
	mpz_sub_ui(u, N, 1);

	gmp_randstate_t state;
	bool flg;

	gmp_randinit_default(state); //srand(time(NULL));

	while (mpz_fdiv_ui(u, 2) == 0)
	{
		mpz_tdiv_q_2exp(u, u, 1); //u >>= 1;
		mpz_add_ui(t, t, 1); //t++
	}

	for (uint32_t cnt = 0; cnt < k; cnt++) {
		mpz_urandomm(a, state, c); //a = 1 + rand() % (c - 1);
		if (!mpz_sgn(a)) {
			mpz_add_ui(a, a, 1);
		}
		mpz_gcd(gcd, a, c);
		if (mpz_cmp_ui(gcd, 1))
			return false;
		mpz_powm(a, a, u, N); //a = (uint64_t)pow((double)a, (double)u) % N;
		if (!mpz_cmp_ui(a, 1)) //if (a == 1)
			continue;
		flg = false;
		for (mpz_set_ui(i, 0); mpz_cmp(i,t) < 0; mpz_add_ui(i, i ,1)) { //for (int i = 0; i < t; i++)
			mpz_mod(a_mod_N, a, N);
			mpz_sub_ui(N_minus_1, N, 1);
			if (!mpz_cmp(a_mod_N, N_minus_1)) { //if (a % N == N - 1) {
				flg = true;
				continue;
			}
			mpz_powm_ui(a, a, 2, N); //a = (a*a) % N;
		}
		if (flg == true)
			continue;
		else
			return false;
	}
	return true;
}
