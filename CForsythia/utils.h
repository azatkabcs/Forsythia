#pragma once
#include <iostream>
using namespace std;

#define Type int

struct out_egcd_t {
	Type d;
	Type u;
	Type v;
};

typedef struct out_egcd_t out_egcd;

Type GCD(Type a, Type b);
void EGCD(Type a, Type b, out_egcd &out);
void dump_egcd(out_egcd);
Type modinv(Type a, Type m);
