/*     Residual Smoothing */

#include"TRI.h"

void smoothing(int N, double ep, vector<vector<double> > &rcs, vector<vector<double> > &rus, vector<vector<double> > &rvs)
{
    vector<vector<double> > Dummyu (N,vector<double>(N, 0));
    vector<vector<double> > Dummyv (N,vector<double>(N, 0));
    vector<vector<double> > Dummyc (N,vector<double>(N, 0));
    
    vector<double> a(N);
    vector<double> b(N);
    vector<double> c(N);
    vector<double> du(N);
    vector<double> dv(N);
    vector<double> dc(N);


    for (int i = 0; i < N; ++i)
    {
        a[i]=0;
        b[i]=0;
        c[i]=0;
        dc[i]=0;
        du[i]=0;
        dv[i]=0;
        for (int j = 0; j < N; ++j)
        {
            Dummyc[i][j]=0;
            Dummyu[i][j]=0;
            Dummyv[i][j]=0;            
        }
    }
   
    a[0]=0;
    b[0]=1+2*ep;
    c[0]=-ep;
    a[N-1]=-ep;
    b[N-1]=1+2*ep;
    c[N-1]=0;
    for (int i = 0; i < N; ++i)
    {
        dc[0]=rcs[i][0];
        du[0]=rus[i][0];
        dv[0]=rvs[i][0];
        for (int j = 1; j < N-1; ++j)
        {
            a[j]=-ep;
            b[j]=1+2*ep;
            c[j]=-ep;
            dc[j]=rcs[i][j];
            du[j]=rus[i][j];
            dv[j]=rvs[i][j];
        }
        dc[N-1]=rcs[i][N-1];
        du[N-1]=rus[i][N-1];
        dv[N-1]=rvs[i][N-1];
        TRI(1,N-1,a,b,c,dc);
        TRI(1,N-1,a,b,c,du);
        TRI(1,N-1,a,b,c,dv);
        dv[0]=rvs[i][0];
        for (int j = 0; j < N; ++j)
        {
            Dummyc[i][j]=dc[j];
            Dummyu[i][j]=du[j];
            Dummyv[i][j]=dv[j];
        }
    }
        
    for (int i = 0; i < N; ++i)
    {
        dc[0]=Dummyc[i][0];
        du[0]=Dummyu[i][0];
        dv[0]=Dummyv[i][0];
        for (int j = 1; j < N-1; ++j)
        {
            a[j]=-ep;
            b[j]=1+2*ep;
            c[j]=-ep;
            dc[j]=Dummyc[i][j];
            du[j]=Dummyu[i][j];
            dv[j]=Dummyv[i][j];
        }
        dc[N-1]=Dummyc[i][N-1];
        du[N-1]=Dummyc[i][N-1];
        dv[N-1]=Dummyc[i][N-1];
        TRI(1,N-1,a,b,c,dc);
        TRI(1,N-1,a,b,c,du);
        TRI(1,N-1,a,b,c,dv);
        for (int j = 0; j < N; ++j)
        {
            rcs[i][j]=dc[j];
            rus[i][j]=du[j];
            rvs[i][j]=dv[j];
        }
    } 
}
