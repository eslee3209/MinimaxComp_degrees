# MinimaxComp_degrees
This algorithm finds optimized degrees for comparison/max/ReLU algorithms using minimax composite polynomial on the RNS-CKKS scheme, which was proposed in https://ieeexplore.ieee.org/document/9517029 and https://eprint.iacr.org/2021/1215. 

## Build
This code requires cmake, NTL, and gmp library. Libraries should be in ```/usr/local/lib```, and header files related to libraries should be in ```/usr/local/include``` or ```/usr/local/include/NTL```. You can perform this algorithm by the following commands in the main directory "MinimaxComp_degrees":
+ ```cmake -S . -B build```
+ ```cd build```
+ ```make```
+ ```./degrees```

## User Setting Parameters
You should set parameters to use in "run_compare.cpp". 
+ alpha: precision parameter 
+ max_factor: max function factor. 1 for comparison operation, and larger than 1 for max/ReLU operation
+ maxdeg: maximum degree. 31 or 63
+ is_comp: true for comparison operation, false for max/ReLU operation
+ level: required level consumption for comparison/max/ReLU operation
+ epsilon: precision parameter

Here, epsilon is equal to max_factor * 2^-alpha by default. You can also use other epsilon for comparison operation.

## ComputeMinDep, ComputeMinMultDegs
You should uncomment ```compute_min_multdepth(RR(alpha), epsilon, maxdeg, is_comp);``` to perform ComputeMinDep and ComputeMinMultDegs algorithms. 

## ComputeMinTimeDegs
You should uncomment ```compute_min_multdepth_update(RR(alpha), epsilon, level, maxdeg, is_comp);``` to perform ComputeMinTimeDegs. Level should be larger than or equal to the minimum level consumption.

The obtained degrees from ComputeMinTimeDegs do not always give better performance than those from ComputeMinMultDegs, so the degrees obtained from ComputeMinTimeDegs may be used only when better performance is given as a result of actual execution. 

If you use another library other than the SEAL library to perform comparison/max/ReLU, you may not be able to see the performance improvement from the ComputeMinTimeDegs algorithm. In this case, it is recommended that CCx.txt files be newly created in a new library/computer environment. Then ComputeMinTimeDegs will give the optimized degrees for your environment.
