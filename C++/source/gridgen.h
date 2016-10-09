#include<math.h>

using namespace std;


// Intializing the grid
void gridgen(int N, double st, double** x, double** y)
{
	double x1[N],y1[N],a1;
	if (st==1)
	{
	    double dx=2*M_PI/(N-1);
	    x1[0]=0;
	    y1[0]=0;
	    for (int i = 1; i < N; ++i)
	    {
	    	x1[i]=x1[i-1]+dx;
	    	y1[i]=y1[i-1]+dx;
	    }
	 }
	else
	{    
	    a1=(M_PI)*(pow((st),(-(N-1)/2)));
	    x1[0]=a1;
	    y1[0]=a1;
	
		    
		 for(int w=1;w<=((N-1)/2); w++)
		 { 	
			x1[w]=a1*(pow((st),w));
			y1[w]=a1*(pow((st),w));
		  }	
		  for(int w=((N+1)/2); w<N; w++)
		  {	
			x1[w]=(2*M_PI)-x1[N-w];
			y1[w]=(2*M_PI)-x1[N-w];
		   }
		    x1[N-1]=2*M_PI;
		    y1[N-1]=2*M_PI;
	    }
	for (int i=0; i<N; i++)
	{
	    for (int j=0; j<N; j++)
	    {
		x[i][j]=x1[j];
		y[i][j]=y1[i];
	    }
	}
}
