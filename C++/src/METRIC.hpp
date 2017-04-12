/*

	Class to do metric transformation. currently it supports 2D transformation. 3D transformation abilities to be added in future
	Created on: April 5, 2017
    Author: Sagar Bhatt

*/


#include <vector>

class METRIC
{
	
public:
	METRIC(){}
	METRIC(int& N, vector<double> &x, vector<double> &y, vector<vector<double> > &zx, vector<vector<double> > &ey, vector<vector<double> > &J);
	
};

METRIC::METRIC(int& N, vector<double> &x, vector<double> &y, vector<vector<double> > &zx, vector<vector<double> > &ey, vector<vector<double> > &J)
{
	vector<vector<double> > G (N, vector<double>(N,  0));
	vector<double> xz (N, 0);
	vector<double> ye (N, 0);  		//note that xe and yz are 0. Hence not computed.
	for (int i=1; i<N-1; i++)
	{
		xz[i]=0.5*(x[i+1]-x[i-1]);
		ye[i]=0.5*(y[i+1]-y[i-1]);
	    
	}
    xz[0]=(x[1]-x[0]);
    ye[0]=(y[1]-y[0]);
    xz[N-1]=(x[N-1]-x[N-2]);
    ye[N-1]=(y[N-1]-y[N-2]);
    

	for (int i=0; i<N; i++)
	{
	    for (int j=0; j<N; j++)
		{
			G[i][j]=xz[i]*ye[j];
			J[i][j]=1/G[i][j];
			zx[i][j] = 1/xz[i];
			ey[i][j] = 1/ye[j];
    	}
	}


}