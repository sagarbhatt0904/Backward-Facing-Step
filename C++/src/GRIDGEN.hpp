/*

	Class to initialize grid points. Currently it supports cartesian coordinate system only.
	Created on: April 5, 2017
    Author: Sagar Bhatt

*/



#include <vector>

class GRIDGEN
{
	double dim_step;	// spatial step size

public:
	GRIDGEN(){}
	GRIDGEN(int &N, double& dim, std::vector<double> &v)
	{
		dim_step=dim/(N-1);
		v.resize(N);
		v[0]=0;
		for (int i = 1; i < N; ++i)
        {
            v[i]=v[i-1]+dim_step;
        }
	}
	
};

/*

write a method for stretched grid

*/