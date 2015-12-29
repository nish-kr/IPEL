//IPEZlwe 
//Inner Product Encryption over Z using lwe for short integer vectors

#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>

using namespace NTL;

struct MPK{
	int n, l;
	ZZ P, V, K;
	ZZ q, m;
	double alp;
	mat_ZZ_p A, U;
};

struct MSK{
	int n, l;
	ZZ P, V, K;
	ZZ q, m;
	double alp;
	mat_ZZ Z;
};

struct SK{
	vec_ZZ sk;
};

struct CT{
	vec_ZZ_p c0;
	vec_ZZ_p c1;
};

void setup(int n, int l, const ZZ& p, const ZZ& v, MPK& mpk, MSK& msk);
SK keygen(const MSK& msk, const ZZ& x);
CT encrypt(const MPK& mpk, const ZZ& y);
ZZ decrypt(const MPK& mpk, const ZZ& x, const SK& sk, const CT& ct);