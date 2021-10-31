#include "func.h"

RR GetInvApproxError(int d, RR t, vector<RR> &X, vector<RR> &Y, long num) {
	if(t<=RR(0) || t>=RR(1)) throw std::out_of_range("t should be in (0,1)");
	return  exptoreal(invexpmaxerr(d, realtoexp(t), X, Y, num));
}
int dep(int deg) // only odd degrees
{
	int d[32]={0,2,3,3,4,4,4,4, 5,5,5,5,5,5,5,5, 6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};	// deg: 1,3,5,7,...,63
	return d[(deg-1)/2];
}
int mult(int deg) // only odd degrees
{
	int m[32]={0,2,3,5,5,6,7,8,8,8,9,9,10,10,11,12, 11,11,11,11,12,12,13,13,14,14,14,14,15,15,16,17};	// deg: 1,3,5,7,...,63
	return m[(deg-1)/2];
}
RR exptoreal(RR x) {
	if(x<0) return pow(static_cast<RR>(2.0),x);
	else return 1.0 - pow(static_cast<RR>(2.0),-x);
}
RR realtoexp(RR x) {
	if(x<=0) {
		cout << "error occur" << endl;
		cout << x << endl;
	}
	if(x>=1) return RR(10000);
	else if(x < 0.5) return log(x)/log(static_cast<RR>(2));
	else return -log(static_cast<RR>(1)-x)/log(static_cast<RR>(2));
	
}
RR invexpmaxerr(long deg, RR expy, vector<RR> &X, vector<RR> &Y, long num) {
	RR expx;
	
	for(long i=1; i<num; i++) {
		if(expy <= Y[i]) {
			// expx ~~ 1.0  or -1.0 < expy < 1.0
			if(X[i-1]<0 && X[i]>0) expx = X[i];
			else if(Y[i-1]<0 && Y[i]>0) expx = X[i];

			// normal case
			else expx = X[i-1] + (X[i]-X[i-1]) * (expy - Y[i-1]) / (Y[i] - Y[i-1]);
			break;
		}

		// expy > 32.0
		if(i==num-1) {
			expx = X[i];
			break;
		}
	}
	return expx;
}
