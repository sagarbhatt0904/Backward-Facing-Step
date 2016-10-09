#include<math.h>

using namespace std;
/*Initializing U,V, Pressure*/
void init(int N,double** x,double** y,double** xvel,double** xvel1,double** yvel,double** yvel1,double** Press,double** dummyu)
{
	

	//function[xvel,yvel,Press,dummyu]=init(N,y)
	/* Initializing U,V, P*/
	int kl=0;
	double xd[N][N];
	for (int i=0; i<N; i++)
	{
	    for (int j=0; j<N; j++)
	    {
	    	xd[i][j]=0;
			xvel1[i][j]=0;
			yvel[i][j]=0;
			Press[i][j]=0;
		}
	}
	for (int j=((N+1)/2)+1; j<N ; j++)
	{
	    xd[1][j]=1-pow((y[j][1]-1),2);
	} 
	for (int j=0; j<(N+1)/2; j++)
	{
		xd[1][j]= xd[1][(N)-kl];
	    kl=kl+1;
	}   

	xd[1][(N+1)/2]=1;
	for(int j=(N+1)/2; j<N; j++)
	{    
	    xvel1[1][j]=xd[1][2*j-N];
	}
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			xvel[i][j]=xvel1[j][i];
			dummyu[i][j]=xvel[i][j];
		}
	}


}