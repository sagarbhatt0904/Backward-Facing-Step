
// Metric calculation

void metric(int N, vector<vector<double> > &x, vector<vector<double> > &y, vector<vector<double> > &zx,  vector<vector<double> > &zy, vector<vector<double> > &ex,  vector<vector<double> > &ey, vector<vector<double> > &JC)
{

vector<vector<double> > G (N, vector<double>(N,  0));
vector<vector<double> > xz (N, vector<double>(N,  0));
vector<vector<double> > ye (N, vector<double>(N,  0));
vector<vector<double> > yz (N, vector<double>(N,  0));
vector<vector<double> > xe (N, vector<double>(N,  0));

for (int i=1; i<N-1; i++)
	{
	    for (int j=0; j<N; j++)
    	    {
		xe[i][j]=0.5*(x[i+1][j]-x[i-1][j]);
		ye[i][j]=0.5*(y[i+1][j]-y[i-1][j]);
	    }
	}
for (int j=0; j<N; j++)
    {
	    xe[0][j]=(x[1][j]-x[0][j]);
	    ye[0][j]=(y[1][j]-y[0][j]);
	}

for (int j=0; j<N; j++)
   {	
	    xe[N-1][j]=(x[N-1][j]-x[N-2][j]);
	    ye[N-1][j]=(y[N-1][j]-y[N-2][j]);
    }

for (int i=0; i<N; i++)
	{
	    for (int j=1; j<N-1; j++)
    	    {
		yz[i][j]=0.5*(y[i][j+1]-y[i][j-1]);
		xz[i][j]=0.5*(x[i][j+1]-x[i][j-1]);
	    }
	 }


for (int i=0; i<N; i++)
	{
	    yz[i][1]=(y[i][2]-y[i][1]);
	    xz[i][1]=(x[i][2]-x[i][1]);
	}
for (int i=0; i<N; i++)
	{
	    yz[i][N-1]=(y[i][N-1]-y[i][N-2]);
	    xz[i][N-1]=(x[i][N-1]-x[i][N-2]);
	}


for (int i=0; i<N; i++)
	{
	    for (int j=1; j<N-1; j++)
    	    {
		G[i][j]=xz[i][j]*ye[i][j]-xe[i][j]*yz[i][j];
		JC[i][j]=1/G[i][j];
    	     }
	}

for (int i=0; i<N; i++)
{
	for (int j=1; j<N-1; j++)
	{
	zx[i][j] = ye[i][j]/G[i][j];
	ex[i][j] = -yz[i][j]/G[i][j];
	zy[i][j] = -xe[i][j]/G[i][j];
	ey[i][j] = xz[i][j]/G[i][j];
	}	
}
for (int i=0; i<N; i++)
{
		JC[i][0]=JC[i][1];
		JC[0][i]=JC[1][i];
		JC[N-1][i]=JC[N-2][i];
		JC[i][N-1]=JC[i][N-2];
}
JC[0][0]=JC[1][1];
}