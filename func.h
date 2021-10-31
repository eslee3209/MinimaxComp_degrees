#pragma once

#include <iostream>
#include <NTL/RR.h>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

using namespace std;
using namespace NTL;

RR GetInvApproxError(int d, RR t, vector<RR> &X, vector<RR> &Y, long num);
int dep(int deg);
int mult(int deg);
RR exptoreal(RR x);
RR realtoexp(RR x);
RR expmaxerr(long deg, RR expx);
RR invexpmaxerr(long deg, RR expy, vector<RR> &X, vector<RR> &Y, long num);