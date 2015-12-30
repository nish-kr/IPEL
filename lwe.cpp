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

void sampleZDiscGaussian(ZZ& ans, RR std_dev){
	//Using plain rejection sampling with support from NTL to compute transcedantal functions
	ZZ maxx = conv<ZZ>(ceil(12*std_dev));
	ZZ temp = 2*maxx+1;
	RR temp1 = inv(2*sqr(std_dev)); // 1/(2*sigma*sigma)
	while(true){
		//Current random value from [-maxx,maxx] intersect Z
		RandomBnd(ans, temp);
		sub(ans,ans,maxx);
		//Accept or reject with apt prob
		RR samp = random_RR();
		RR temp2 = exp(-1*(conv<RR>(sqr(ans)))*temp1);
		if (samp < temp2)
			return ;
	}
}

void gau_sampleZ(vec_ZZ& ans, int dim, RR stddev){
	ans.SetLength(dim);
	for(int i=0;i<dim;i++){
		sampleZDiscGaussian(ans[i], stddev);
	}
}

ZZ findMini(const MPK& mpk, const ZZ_p& muPr){
	ZZ i(0), j = mpk.K-1;
	ZZ target = conv<ZZ>(muPr);
	ZZ di = mpk.q/mpk.K;
	while(i<j){
		ZZ mid = (i+j)/2;
		ZZ midp = mid+1;
		ZZ dist = abs(di*mid-target);
		ZZ distp = abs(di*midp-target);
		if (4*dist <= di){
			return mid;
		}
		else if (4*distp <= di){
			return midp;
		}
		else if (distp < dist){
			i = mid;
		}
		else if (distp > dist){
			j = midp;
		}
	}
	// ?? i>j ?
	return i;
}

void setup(int n, int l, const ZZ& p, const ZZ& v, MPK& mpk, MSK& msk, int logq){
	//remember to initialize ZZ_p in setup
	mpk.P = msk.P = p;
	mpk.V = msk.V = v;
	mpk.K = msk.K = l*p*v;
	
	mpk.n = msk.n = n;
	mpk.l = msk.l = l;
	mpk.m = msk.m = 4*n*logq;

	int temp1 = ceil(pow((double)n,1.5));
	temp1 *= logq;
	ZZ qlow = temp1*sqr(mpk.K);
	ZZ ql = NextPrime(qlow);
	mpk.q = msk.q = ql;

	//Initializing ZZ_p with q
	ZZ_p::init(ql);

	ll temp2 = (ll)n*n*n*n*(mpk.m)*((ll)mpk.m);
	double temp3 = pow(log2((double)n),3.5);
	double alp = 1.0/(temp2*temp3);
	mpk.alp = msk.alp = alp;


	//compute sigma1
	double temp4 = pow(((double)n)*log2((double)mpk.m), 0.5);
	int temp5 = ceil(pow((double)mpk.m,0.5));
	RR sigma1 = conv<RR>(temp4)*conv<RR>(max(conv<ZZ>(temp5),mpk.K));

	//compute sigma2
	double temp6 = pow((double)n,3.5)*pow((double)mpk.m,0.5)*pow(log2((double)mpk.m),2.5);
	ZZ temp7 = max(conv<ZZ>(mpk.m), sqr(mpk.K));
	RR sigma2 = conv<RR>(temp6)*conv<RR>(temp7);

	//fill Z
	msk.Z.SetDims(l,mpk.m);
	for(int i=0;i<l;i++){
		for(int j=0;j<mpk.m;j++){
			if (j<(mpk.m)/2){
				sampleZDiscGaussian(msk.Z[i][j], sigma1);
			}
			else if (j == i+(mpk.m)/2){
				sampleZDiscGaussian(msk.Z[i][j], sigma2);
				//This was according to centre 0. Now shift to centre to 1.
				//???
				msk.Z[i][j] += 1;
			}
			else{
				sampleZDiscGaussian(msk.Z[i][j], sigma2);
			}
		}
	}

	//fill A
	mpk.A.SetDims(mpk.m,n);
	for(int i=0;i<mpk.m;i++){
		for(int j=0;j<n;j++){
			random(mpk.A[i][j]);
		}
	}

	//fill U
	mul(mpk.U, conv<mat_ZZ_p>(msk.Z), mpk.A);
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

	//NOTE: converting q to double .. ??
	vec_ZZ e0,e1;
	gau_sampleZ(e0, mpk.m, conv<RR>(mpk.alp)*conv<RR>(mpk.q));
	gau_sampleZ(e1, mpk.l, conv<RR>(mpk.alp)*conv<RR>(mpk.q));

	vec_ZZ_p s;
	s.SetLength(mpk.n);
	for(int i=0; i<mpk.n; i++){
		random(s[i]);
	}

	mul(ans.c0, mpk.A, s);
	ans.c0 += conv<vec_ZZ_p>(e0);

	mul(ans.c1, mpk.U, s);
	ans.c1 += conv<vec_ZZ_p>(e1);
	ans.c1 += (conv<ZZ_p>(mpk.q/mpk.K))*conv<vec_ZZ_p>(y);
}

ZZ decrypt(const MPK& mpk, const vec_ZZ& x, const SK& sk, const CT& ct){
	assert(x.length() == mpk.l);
	//assert each element of x is within V.
	doassert(x, mpk.V);

	ZZ_p temp1, temp2;
	InnerProduct(temp1, ct.c1, conv<vec_ZZ_p>(x));
	InnerProduct(temp2, ct.c0, conv<vec_ZZ_p>(sk.sk));
	sub(temp1,temp1,temp2);
	ZZ ans = findMini(mpk, temp1);
	return ans;
}

int main(){
	ZZ x(29);
	ZZ temp(5);
	ZZ_p::init(temp);
	ZZ_p y = conv<ZZ_p>(x);
	cout<<y<<endl;
}