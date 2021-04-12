/* GCD, Extended GCD and Modular inversion implementation */
#include "utils.h"

/* grand common divisor */
Type GCD(Type a, Type b) {
	Type r;

	while (!b == 0) {
		r = a % b;
		a = b;
		b = r;
	}
	return a;
}

/* extended grand common divisor */
void EGCD(Type a, Type b, out_egcd &out) {
	Type v1, v3, q, t1, t3;
	out.u = 1;
	out.d = a;

	if (b == 0) {
		out.v = 0;
#ifdef DEBUG
		dump_egcd(out);
#endif
		return;
	}
	else {
		v1 = 0;
		v3 = b;
	}

	while (v3 != 0) {
		q = out.d / v3;
		t3 = out.d % v3;
		t1 = out.u - q * v1;
		out.u = v1;
		out.d = v3;
		v1 = t1;
		v3 = t3;
	}
	out.v = (out.d - a * out.u) / b;
#ifdef DEBUG
	dump_egcd(out);
#endif
}

void dump_egcd(out_egcd out) {
	cout << "Extended GCD" << endl;
	cout << "d = " << out.d << " " << "u = " << out.u << " " << "v = " << out.v << endl;
}

Type modinv(Type a, Type m) {
	out_egcd out = { 0 };
	EGCD(a % m, m, out);

	if (out.d != 1)
		cout << "Error! Modular inverse does not exist!" << endl;
	else
		return out.u + m;

	return -1;
}
