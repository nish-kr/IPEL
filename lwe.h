//IPEZlwe 
//Inner Product Encryption over Z using lwe for short integer vectors

#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/RR.h>

using namespace NTL;

typedef long long ll;

struct MPK{
	int n, l, m; // NOTE: taking m as an int
	ZZ P, V, K;
	ZZ q;
	double alp; // NOTE: Not sure .. double or RR
	mat_ZZ_p A, U;
};

struct MSK{
	int n, l, m; // NOTE: taking m as an int
	ZZ P, V, K;
	ZZ q;
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

void setup(int n, int l, const ZZ& p, const ZZ& v, MPK& mpk, MSK& msk, int logq);
SK keygen(const MSK& msk, const vec_ZZ& x);
CT encrypt(const MPK& mpk, const vec_ZZ& y);
ZZ decrypt(const MPK& mpk, const vec_ZZ& x, const SK& sk, const CT& ct);