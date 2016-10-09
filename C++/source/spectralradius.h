#include<math.h>

using namespace std;

//Computing spectral radius
void spectralradius(int N,double** JC,double** U,double** V,double** g11,double** g22, double** rho1, double** rho2)
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
