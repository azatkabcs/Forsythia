// CForsythia.cpp: определяет точку входа для приложения.
//

#include "CForsythia.h"

using namespace std;
vector<uint32_t> smallPrimes;

int main()
{
	mpz_t x, y, eA, eB, f;
	mpz_inits(x, y, p_global, eA, eB, f, NULL);
	mpz_ui_pow_ui(eA, 2, 137);
	mpz_ui_pow_ui(eB, 3, 84);
	mpz_set_ui(f, 1);
	mpz_mul(p_global, eA, eB);
	mpz_mul(p_global, p_global, f);
	mpz_sub_ui(p_global, p_global, 1);
	mpz_set_str(x, "1697515660496260333664277950769004218366836292978277269037599020171750501119102405", 10);
	GFp2Element A(x, y);
	GFp2Element B(y, x);

	gmp_printf("global p %Zd\n", p_global);
	if (A != B)
		cout << "here" << endl;
	cout << A << endl;
	gmp_printf("%Zd\n", x);
	return 0;
}
