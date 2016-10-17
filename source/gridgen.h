#include<math.h>


// Intializing the grid
void gridgen(int N, double st, double xmax, double ymax, double** x, double** y)
{
    double x1[N],y1[N],a1;
	if (st==1)
	{
	    double dx=xmax/(N-1);
        double dy=ymax/(N-1);
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
        double dx=(xmax*(1-st)/(1-(pow(st,(N-1)))));
		double dy=(ymax*(1-st)/(1-(pow(st,(N-1)))));
        x1[0]=0;
        y1[0]=0;
        for (int i = 1; i < N; ++i)
	    {
            x1[i]=x1[i-1]+(pow(st,(i-2)))*dx;
            y1[i]=y1[i-1]+(pow(st,(i-2)))*dy;
        }
        x1[N-1]=xmax;
        y1[N-1]=ymax;
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