#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

using namespace NTL;
using namespace std;

int main(){
	mat_ZZ x;
	int m=2,n=3;
	x.SetDims(m,n);
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			x[i][j] = m*i+n*j;
			cout<<x[i][j]<<" ";
		}
		cout<<endl;
	}
}