#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_ZZ.h>

using namespace NTL;

struct MSK{
	vec_ZZ_p s;
	ZZ p;
	ZZ_p g;
	int l; 		//the vector length being used
	int sec; 	//the security parameter
};

struct MPK{
	vec_ZZ_p h;
	ZZ p;
	ZZ_p g;
	int l; 		//the vector length being used
	int sec; 	//the security parameter
};

struct SK{
	ZZ sk;
};

struct CT{
	ZZ_p c0;
	vec_ZZ_p c1;
};


void setup(MPK& mpk, MSK& msk, int sec, int l, int pbitl);//pbitl represents the length of prime to use - it defaults to 512
void encrypt(CT& ct, const MPK& mpk, const vec_ZZ_p msg); //msg must have exactly l elements
void keyder(SK& seckey, const MSK& msk, const vec_ZZ_p& y);
ZZ_p decrypt(const MPK& mpk, const CT& ct, SK& seckey, const vec_ZZ_p& y);