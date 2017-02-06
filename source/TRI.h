

// Tridiagonal Matrix Solver
void TRI(int ibeg,int iend,vector<double> &aa,vector<double> &bb,vector<double> &cc,vector<double> &dd)
{
	for(int i=ibeg; i<iend; i++)
	{
	    double r=aa[i]/bb[i-1];
	    bb[i]=bb[i]-r*cc[i-1];
	    dd[i]=dd[i]-r*dd[i-1];  

	}
	dd[iend]=dd[iend]/bb[iend];

	for(int i=iend-2;i>=ibeg;i--)
	{
	    dd[i]=(dd[i]-cc[i]*dd[i+1]/bb[i]);
	}
}
