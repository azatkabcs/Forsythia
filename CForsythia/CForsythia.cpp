// CForsythia.cpp: определяет точку входа для приложения.
//

#include "CForsythia.h"

using namespace std;
vector<uint32_t> smallPrimes;

int main()
{
	mpz_t x, y, eA, eB, f, index;
	mpz_inits(x, y, p_global, eA, eB, f, index, NULL);
	mpz_ui_pow_ui(eA, 2, 137);
	mpz_ui_pow_ui(eB, 3, 84);
	mpz_set_ui(f, 1);
	mpz_mul(p_global, eA, eB);
	mpz_mul(p_global, p_global, f);
	mpz_sub_ui(p_global, p_global, 1);
	mpz_set_str(x, "1697515660496260333664277950769004218366836292978277269037599020171750501119102405", 10);
	mpz_set_ui(index, 1000);
	
	GFp2Element A(x, y);
	GFp2Element B(y, x);

	MontgomeryCurve E0(A, f);
	MontgomeryPoint point(A, A, A, f);

	cout << point << endl;
	cout << point.mul2() << endl;
	cout << point.mul3() << endl;
	cout << point.mul2e(index) << endl;
	cout << point.mul3e(index) << endl;
	mpz_clears(x, y, p_global, eA, eB, f, index, NULL);
	return 0;
}
