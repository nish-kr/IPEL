#include "ddh.h"
#include <iostream>
#include <cassert>
#include <cstdlib>

using namespace NTL;

ZZ_p find_primitive_root(ZZ p, ZZ q){
	//p-1 = 2q
	ZZ_p temp1, temp2;
	while(true){
		random(temp1);
		if (temp1 == 0) continue;
		power(temp2, temp1, 2);
		if (temp2 != 1){
			power(temp2, temp1, q);
			if (temp2 != 1){
				break;
			}
		}
	}
	return temp1;
}


void compute_disc_log(ZZ_p& ans, const MPK& mpk, const ZZ_p& x){
	ans = 0;
	if (x==1) return;
	ZZ_p temp = mpk.g;
	ans += 1;
	while(temp != x && ans != 0){
		ans += 1;
		temp *= mpk.g;
		// std::cout<<ans<<" "<<temp<<std::endl;
	}
}


void setup(MPK& mpk, MSK& msk, int sec, int l, int pbitl){
	std::cout<<"Executing Setup--\n";
	ZZ p, q;
	GenGermainPrime(q, pbitl, sec);
	p = 2*q+1;
	ZZ_p::init(p);
	std::cout<<"Executing find_primitive_root--\n";
	ZZ_p g = find_primitive_root(p, q);
	msk.s.SetLength(l);
	mpk.h.SetLength(l);
	for(int i=0;i<l;i++){
		random(msk.s[i]);
		power(mpk.h[i], g, conv<ZZ>(msk.s[i]));
	}
	msk.p = mpk.p = p;
	msk.g = mpk.g = g;
	msk.l = mpk.l = l;
	msk.sec = mpk.sec = sec;
	std::cout<<"Finished Setup--\n";
}


void encrypt(CT& ct, const MPK& mpk, const vec_ZZ_p msg){
	assert(msg.length() == mpk.l);
	std::cout<<"Executing encrypt--\n";
	ZZ_p r_zpp;
	random(r_zpp);
	ZZ r_z = conv<ZZ>(r_zpp);
	std::cout<<"Encrypt r chosen = "<<r_zpp<<"\n";
	ct.c1.SetLength(mpk.l);
	ZZ_p temp1, temp2;
	for(int i=0,len=mpk.l;i<len;i++){
		power(temp1, mpk.h[i], r_z);
		power(temp2, mpk.g, conv<ZZ>(msg[i]));
		mul(ct.c1[i], temp1, temp2);
		std::cout<<"ct.c1["<<i<<"] = "<<ct.c1[i]<<"\n";
	}
	power(ct.c0, mpk.g, r_z);
	std::cout<<"ct.co = "<<ct.c0<<"\n";
	std::cout<<"Finished encrypt--\n";
}


void keyder(SK& seckey, const MSK& msk, const vec_ZZ_p& y){
	std::cout<<"Executing keyder--\n";
	InnerProduct(seckey.sk, conv<vec_ZZ>(msk.s), conv<vec_ZZ>(y));
	std::cout<<"Secret key = "<<seckey.sk<<std::endl;
}


ZZ_p decrypt(const MPK& mpk, const CT& ct, SK& seckey, const vec_ZZ_p& y){
	std::cout<<"Executing decrypt--\n";
	ZZ_p temp1 = ZZ_p(1), temp2, temp3;
	for(int i=0,len=mpk.l;i<len;i++){
		power(temp2, ct.c1[i], conv<ZZ>(y[i]));
		temp1 *= temp2;
	}
	std::cout<<"temp1 = "<<temp1<<std::endl;
	power(temp2, ct.c0, (seckey.sk));
	std::cout<<"temp2 = "<<temp2<<std::endl;

	div(temp3, temp1, temp2);
	std::cout<<"temp3 = "<<temp3<<std::endl;

	std::cout<<"Executing compute_disc_log--\n";
	compute_disc_log(temp1, mpk, temp3);
	std::cout<<"Finished decrypt--\n";
	return temp1;
}


void printMPK(const MPK& mpk){
	std::cout<<"--------------------\n";
	std::cout<<"Printing MPK == \n";
	std::cout<<"MPK.p = "<<mpk.p<<"\n";
	std::cout<<"MPK.g = "<<mpk.g<<"\n";
	std::cout<<"MPK.l = "<<mpk.l<<"\n";
	std::cout<<"MPK.sec = "<<mpk.sec<<"\n";
	for(int i=0,len=mpk.l;i<len;i++){
		std::cout<<"MPK.h["<<i<<"] = "<<mpk.h[i]<<"\n";
	}
	std::cout<<"--------------------\n";
}


void printMSK(const MSK& msk){
	std::cout<<"--------------------\n";
	std::cout<<"Printing MSK == \n";
	std::cout<<"MSK.p = "<<msk.p<<"\n";
	std::cout<<"MSK.g = "<<msk.g<<"\n";
	std::cout<<"MSK.l = "<<msk.l<<"\n";
	std::cout<<"MSK.sec = "<<msk.sec<<"\n";
	for(int i=0,len=msk.l;i<len;i++){
		std::cout<<"MSK.s["<<i<<"] = "<<msk.s[i]<<"\n";
	}
	std::cout<<"--------------------\n";
}


int main(){
	MPK mpk;
	MSK msk;
	int l = 2;
	setup(mpk,msk,80,l,512);
	printMSK(msk);
	printMPK(mpk);
	vec_ZZ_p x,y;
	x.SetLength(l);
	y.SetLength(l);
	ZZ bans(0);
	for(int i=0;i<l;i++){
		// random(x[i]);
		// random(y[i]);
		x[i] = rand()%100;
		y[i] = rand()%100;
		std::cout<<"x["<<i<<"]="<<x[i]<<", y["<<i<<"]="<<y[i]<<"\n";
		bans += (conv<ZZ>(x[i])*conv<ZZ>(y[i]));
		// bans %= conv<ZZ>(mpk.p-1);
	}
	std::cout<<"bans = "<<bans<<std::endl;
	CT ct;
	SK sk;
	encrypt(ct,mpk,x);
	keyder(sk,msk,y);
	ZZ_p ans = decrypt(mpk,ct,sk,y);
	bans %= (mpk.p-1);
	std::cout<<ans<<" "<<bans<<"\n"<<(ans==conv<ZZ_p>(bans))<<std::endl;
}