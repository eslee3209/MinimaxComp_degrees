#include "optimized_degrees.h"

// IEEE Access version
void compute_min_multdepth_update(RR alpha, RR epsilon, int depth, long maxdeg, bool is_comp)
{
	// variables
	string filename;
	RR tau;
	RR target;
	vector<int> lt;
	vector<int>::iterator iter;
	RR temp;
	size_t mindep, mintime;
	size_t t_max=300, n_max=30;	// maxdeg should be odd!
	size_t first_level = depth;
	if(is_comp == false) depth -= 1;

	vector<vector<RR>> u(t_max+1, vector<RR>(n_max+1, RR(0)));
	vector<vector<vector<int>>> V(t_max+1, vector<vector<int>>(n_max+1, vector<int>(0)));
	RR::SetPrecision(300);
	RR::SetOutputPrecision(20);	

	// import StoreME result and runtime text files
	vector<ifstream> E(maxdeg+1), C(maxdeg+1);		// E[0] ~ E[maxdeg]. E of odd index(3~maxdeg) is only used.
	for(size_t deg=3; deg<=maxdeg; deg+=2)
	{
		E[deg].open("../text/E" + to_string(deg) + ".txt");	if(!E[deg].is_open()) throw std::runtime_error("file is not open");
		C[deg].open("../text/CC" + to_string(deg) + ".txt");	if(!C[deg].is_open()) throw std::runtime_error("file is not open");
	}

	// import data from StoreME result and runtime text files
	vector<vector<RR>> X(maxdeg+1, vector<RR>(0)), Y(maxdeg+1, vector<RR>(0));
	vector<vector<vector<size_t>>> time(maxdeg+1, vector<vector<size_t>>(n_max+1, vector<size_t>(n_max+1,0)));
	RR expx;
	vector<size_t> num(maxdeg+1,0);
	for(size_t deg=3; deg<=maxdeg; deg+=2) E[deg] >> num[deg];		// number of pairs
	for(size_t deg=3; deg<=maxdeg; deg+=2)
	{
		RR temp;
		for(size_t i=0; i<num[deg]; i++) 
		{
			E[deg] >> temp;
			X[deg].emplace_back(temp);
			E[deg] >> temp;
			Y[deg].emplace_back(temp);
		}
		for(size_t i=0; i<=n_max; i++)
		{
			for(size_t j=0; j<=n_max; j++)
			{
				C[deg] >> time[deg][i][j];
				time[deg][i][j] = (time[deg][i][j]+50)/100;
			}
		}
	}
	
	// f(m,n,t), G(m,n,t) evaluation
	// parameter setting
	target = (RR(1)-epsilon)/(RR(1)+epsilon);	// the first domain 
	cout << "------------------------------------" << endl;
	cout << "alpha: " << alpha << endl;
	cout << "epsilon: " << epsilon << endl;
	long a=0;
	while(RR(a)+1<alpha+0.5) a++;		// a is an integer that is greater than or equal to alpha

	// evaluates u_tau(m,n), G_tau(m,n)
	for(size_t m=0; m<=t_max; m++) {
		for(size_t n=0; n<=n_max; n++) {
			V[m][n].clear();
//			cout << a << ": " << m << ", " << n << " Round" << endl;

			// m<=1 or n<=1 case
			if(m<=1 || n<=1) {
				u[m][n] = pow(RR(2.0),RR(1)-alpha);
			}
			
			// m>=2 and n>=2 case
			else {
				int j=0;
				RR max = RR(-1);
				for(size_t k=1; 2*k+1<=maxdeg; k++) {
					if(time[2*k+1][first_level][n] <= m && dep(2*k+1) <= n)
					{
						temp = GetInvApproxError(2*k+1,u[m-time[2*k+1][first_level][n]][n-dep(2*k+1)],X[2*k+1],Y[2*k+1],num[2*k+1]);
						
						if(temp > max) {
							j=k;
							max = temp;
						}
					}
				}
				if(max>0)
				{
					lt = V[m-time[2*j+1][first_level][n]][n-dep(2*j+1)];
					u[m][n] = max;
					V[m][n].emplace_back(2*j+1);	
					for(iter=lt.begin(); iter!=lt.end(); iter++) V[m][n].emplace_back(*iter);		
				}
				else u[m][n] = pow(RR(2.0),RR(1)-alpha);
			}
		}
	}

	// ComputeMinTimeDegs
	for(mintime=0; mintime<=t_max; mintime++){
		if(u[mintime][depth] >= target) break;
		if(mintime == t_max) cout << "failure" << endl;
	}
	cout << "mintime: " << mintime << endl;
	if(is_comp == true)	cout << "depth: " << depth << endl;
	else if(is_comp == false) cout << "depth: " << depth+1 << endl;
	cout << 1-u[mintime][depth] << endl;		// This should be larger than epsilon

	lt = V[mintime][depth];
	for(iter=lt.begin(); iter!=lt.end(); iter++) {
		cout << (*iter) << " ";
	}
	cout << endl << endl;
}
// IEEE TDSC version
void compute_min_multdepth(RR alpha, RR epsilon, long maxdeg, bool is_comp){

	// variables
	string filename;
	RR tau;
	RR target;
	vector<int> lt;
	vector<int>::iterator iter;
	RR temp;
	int minmult, mindep;
	// size_t maxdeg=63, m_max=70, n_max=40;	// maxdeg should be odd!
	size_t m_max=70, n_max=40;	// maxdeg should be odd!
	vector<vector<RR>> h(m_max+1, vector<RR>(n_max+1, RR(0)));
	vector<vector<vector<int>>> g(m_max+1, vector<vector<int>>(n_max+1, vector<int>(0)));
	RR::SetPrecision(300);
	RR::SetOutputPrecision(20);	

	// import StoreME result text files
	vector<ifstream> E(maxdeg+1);		// E[0] ~ E[maxdeg]. E of odd index(3~maxdeg) are only used
	for(size_t deg=3; deg<=maxdeg; deg+=2)
	{
		E[deg].open("../text/E" + to_string(deg) + ".txt");
		if(!E[deg].is_open()) throw std::runtime_error("file is not open");
	}

	// import data from StoreME result text files
	vector<vector<RR>> X(maxdeg+1, vector<RR>(0)), Y(maxdeg+1, vector<RR>(0));	// X, Y of odd index(3~maxdeg) are only used
	RR expx;
	vector<size_t> num(maxdeg+1,0);
	for(size_t deg=3; deg<=maxdeg; deg+=2) E[deg] >> num[deg];		// number of pairs is stored in each text file "En.txt"
	for(size_t deg=3; deg<=maxdeg; deg+=2)
	{
		RR temp;
		for(size_t i=0; i<num[deg]; i++) 
		{
			E[deg] >> temp;
			X[deg].emplace_back(temp);
			E[deg] >> temp;
			Y[deg].emplace_back(temp);
		}
	}
	
	// h_tau(m,n), G_tau(m,n) evaluation
	// parameter setting
	target = (static_cast<RR>(1)-epsilon)/(static_cast<RR>(1)+epsilon);	// the first domain 
	cout << "--- parameter ---" << endl;
	cout << "alpha: " << alpha << endl;
	cout << "epsilon: " << epsilon << endl << endl;
	long a=0;
	while(RR(a)+1<alpha+0.5) a++;	// a is an integer that is greater than or equal to alpha, i.e., \lceil alpha \rceil

	// evaluates h_tau(m,n), G_tau(m,n)
	for(size_t m=0; m<=m_max; m++) {
		for(size_t n=0; n<=n_max; n++) {
			g[m][n].clear();
//			cout << a << ": " << m << ", " << n << " Round" << endl;

			// m<=1 or n<=1 case
			if(m<=1 || n<=1) {
				h[m][n] = pow(static_cast<RR>(2.0),static_cast<RR>(1)-alpha);
			}
			
			// m>=2 and n>=2 case
			else {
				int j=0;
				RR max = static_cast<RR>(0);
				for(size_t k=1; mult(2*k+1) <= m && dep(2*k+1) <= n && 2*k+1<=maxdeg; k++) {	
					temp = GetInvApproxError(2*k+1,h[m-mult(2*k+1)][n-dep(2*k+1)],X[2*k+1],Y[2*k+1],num[2*k+1]);
					
					if(temp > max) {
						j=k;
						max = temp;
					}
				}
				lt = g[m-mult(2*j+1)][n-dep(2*j+1)];
				h[m][n] = max;
				g[m][n].emplace_back(2*j+1);	
				for(iter=lt.begin(); iter!=lt.end(); iter++)  g[m][n].emplace_back(*iter);		
			}
		}
	}

	// ComputeMinDep
	cout << "--- ComputMinDep ---" << endl;
	for(mindep=0; mindep<=n_max; mindep++){
		if(h[m_max][mindep] >= target) break;
		if(mindep == n_max) cout << "failure" << endl;
	}
	if(is_comp == true)	cout << "mininum depth: " << mindep << endl << endl;
	else if(is_comp == false) cout << "mininum depth: " << mindep+1 << endl << endl;

	// ComputeMinMultDegs
	cout << "--- ComputMinMultDegs ---" << endl;
	int total_minmult = 1000;
	for(int dep = mindep; dep <= mindep+10; dep++)
	{
		for(minmult=0; minmult<=m_max; minmult++){
			if(h[minmult][dep] >= target) break;
			if(minmult == m_max) cout << "failure" << endl;
		}
		if(minmult == total_minmult) break;
		total_minmult = minmult;

		if(is_comp == true)	cout << "depth: " << dep << endl;
		else if(is_comp == false) cout << "depth: " << dep+1 << endl;
		cout << "minmult: " << minmult << endl;
		cout << 1-h[minmult][dep] << endl;		// This should be larger than epsilon
		cout << "degs: ";
		lt = g[minmult][dep];
		for(iter=lt.begin(); iter!=lt.end(); iter++) {
			cout << (*iter) << " ";
		}
		cout << endl << endl;
	}

}