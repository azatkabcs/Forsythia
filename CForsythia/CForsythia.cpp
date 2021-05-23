#include "CForsythia.h"
#include <chrono>

using namespace std;
vector<uint32_t> smallPrimes;

void init_p_global(parameters params) {
	mpz_t eA, eB, f;
	mpz_inits(eA, eB, f, p_global, NULL);
	mpz_ui_pow_ui(eA, 2, params.eA);
	mpz_ui_pow_ui(eB, 3, params.eB);
	mpz_set_ui(f, params.f);
	mpz_mul(p_global, eA, eB);
	mpz_mul(p_global, p_global, f);
	mpz_sub_ui(p_global, p_global, 1);
	mpz_clears(eA, eB, f, NULL);
}

GFp2Element init_A(parameters params) {
	mpz_t a, b;
	mpz_init_set_str(a, params.A.first, 10);
	mpz_init_set_str(b, params.A.second, 10);
	GFp2Element res(a, b);
	mpz_clears(a, b, NULL);
	return res;
}

GFp2Element init_point(point& x, point& z) {
	mpz_t x1, x2, z1, z2;
	mpz_inits(x1, x2, z1, z2, NULL);
	mpz_set_str(x1, x.first, 10);
	mpz_set_str(x2, x.second, 10);
	mpz_set_str(z1, z.first, 10);
	mpz_set_str(z2, z.second, 10);
	GFp2Element X(x1, x2);
	GFp2Element Z(z1, z2);
	mpz_clears(x1, x2, z1, z2, NULL);
	return GFp2Element(X / Z);
}

void init_pk(mpz_t sk, int base, uint32_t degree) {
#if 1
	mpz_t num;
	mpz_init(num);
	mpz_ui_pow_ui(num, base, degree);
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
	mpz_urandomm(sk, state, num);
	mpz_clear(num);
	gmp_randclear(state);
#else
	if (base == 2) {
		const char num[] = "80057929638325553428237217364670151799135";
		mpz_set_str(sk, num, 10); 
	} else {
		const char num[] = "1908005428125941001067070257899517329963";
		mpz_set_str(sk, num, 10);
	}
#endif
}

void SIDH(parameters params) {
	init_p_global(params);
	GFp2Element A = init_A(params);
	MontgomeryCurve E0(A, GFp2Element(1));
	GFp2Element xp2 = init_point(params.xp2, params.zp2);
	GFp2Element xq2 = init_point(params.xq2, params.zq2);
	GFp2Element yp2 = init_point(params.yp2, params.zp2);
	GFp2Element yq2 = init_point(params.yq2, params.zq2);
	GFp2Element l((yp2 + yq2) / (xp2 - xq2));
	GFp2Element xr2(l * l - (xp2 + xq2) - A);
	GFp2Element xp3 = init_point(params.xp3, params.zp3);
	GFp2Element xq3 = init_point(params.xq3, params.zq3);
	GFp2Element yp3 = init_point(params.yp3, params.zp3);
	GFp2Element yq3 = init_point(params.yq3, params.zq3);

	l = (yp3 + yq3) / (xp3 - xq3);
	GFp2Element xr3(l * l - (xp3 + xq3) - A);
	cout << "A = " << A << endl;
	cout << "E0 = " << E0 << endl;

	mpz_t sk2, eA, sk3, eB;
	mpz_inits(sk2, eA, sk3, eB);
	mpz_set_ui(eA, params.eA);
	mpz_set_ui(eB, params.eB);
	GFp2Element pkA[3], pkB[3], jAlice, jBob;

	auto start_time = std::chrono::steady_clock::now();

	init_pk(sk2, 2, params.eA);
	gmp_printf("%s: %Zd\n", "Alice secret key:", sk2);
	isogen2(E0, sk2, eA, xp2, xq2, xr2, xp3, xq3, xr3, pkA[0], pkA[1], pkA[2]);
	cout << "Alice public key: " << endl;
	cout << pkA[0] << endl << pkA[1] << endl << pkA[2] << endl;

	init_pk(sk3, 3, params.eB);
	gmp_printf("%s: %Zd\n", "Bob secret key:", sk3);
	isogen3(E0, sk3, eB, xp2, xq2, xr2, xp3, xq3, xr3, pkB[0], pkB[1], pkB[2]);
	cout << "Bob public key: " << endl;
	cout << pkB[0] << endl << pkB[1] << endl << pkB[2] << endl;

	jAlice = isoex2(sk2, eA, pkB);
	jBob = isoex3(sk3, eB, pkA);
	cout << "jAlice: " << endl;
	cout << jAlice << endl;
	cout << "jBob: " << endl;
	cout << jBob << endl;

	auto end_time = std::chrono::steady_clock::now();
    auto elapsed_ns = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Elapsed time: "<< elapsed_ns.count() << " ms\n";

	mpz_clears(sk2, eA, sk3, eB, p_global, NULL);
}

int main()
{
	cout << "************************ Set Forsythia80 ************************" << endl;
	SIDH(forsythia80);

	cout << "************************ Set Forsythia128 ************************" << endl;
	SIDH(forsythia128);

	cout << "************************ Set Forsythia256 ************************" << endl;
	SIDH(forsythia256);

	return 0;
}
