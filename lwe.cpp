#include "lwe.h"
#include <iostream>
#include <cassert>
#include <cstdlib>

using namespace std;

void doassert(const vec_ZZ& x, const ZZ& upper){
	int l = x.length();
	for(int i=0; i<l; i++){
		assert(x[i]<upper);
	}
}

void setup(int n, int l, const ZZ& p, const ZZ& v, MPK& mpk, MSK& msk){

}

SK keygen(const MSK& msk, const vec_ZZ& x){
	assert(x.length() == msk.l);
	//assert each element of x is within V.
	doassert(x, msk.V);
	SK ans;
	mul(ans.sk, x, msk.Z);
	return ans;
}

CT encrypt(const MPK& mpk, const vec_ZZ& y){
	assert(y.length() == mpk.l);
	//assert each element of y is within P.
	doassert(y, mpk.P);
	CT ans;
	
}

ZZ decrypt(const MPK& mpk, const vec_ZZ& x, const SK& sk, const CT& ct){

}

int main(){
	ZZ x(29);
	ZZ temp(5);
	ZZ_p::init(temp);
	ZZ_p y = conv<ZZ_p>(x);
	cout<<y<<endl;
}