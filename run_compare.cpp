#include <iostream>
#include <NTL/RR.h>
#include <cmath>
#include "optimized_degrees.h"

using namespace std;
using namespace NTL;

int main() {

	// user setting parameters
	long alpha = 20;		// precision parameter alpha
	long max_factor = 21;	// max_factor = 1 for comparison operation. max_factor > 1 for max/ReLU operation
	long maxdeg = 63;		// 31 or 63
	bool is_comp = false;	// true: comparison operation, false: max/ReLU operation
	long level = 22;		// total level consumption D. just for compute_min_multdepth_update

	// ComputeMinMultDegs or ComputeMinTimeDegs
	// compute_min_multdepth(RR(alpha), max_factor*pow(RR(2),RR(-alpha)), maxdeg, is_comp);
	compute_min_multdepth_update(RR(alpha), max_factor*pow(RR(2),RR(-alpha)), level, maxdeg, is_comp);



	return 0;
}
