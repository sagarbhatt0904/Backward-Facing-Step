#include<math.h>

using namespace std;

//Computing spectral radius
void spectralradius(int N,vector<vector<double> > &JC,vector<vector<double> > &U,vector<vector<double> > &V,vector<vector<double> > &g11,vector<vector<double> > &g22, vector<vector<double> > &rho1, vector<vector<double> > &rho2)
{
	for (int i =0; i<N; i++)
	{
		for (int j =0; j<N; j++)
		{
			rho1[i][j]=(1/JC[i][j])*(abs(U[i][j])+sqrt((pow(U[i][j],2)+g11[i][j])));
			rho2[i][j]=(1/JC[i][j])*(abs(V[i][j])+sqrt((pow(V[i][j],2)+g22[i][j])));
	    }
	}
}
