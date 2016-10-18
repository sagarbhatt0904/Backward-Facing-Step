#include<math.h>
#include"Diss.h"
#include"spectralradius.h"
#include"smoothing.h"


void RHS(int N,vector<vector<double> > &JC,vector<vector<double> > &xvel,vector<vector<double> > &yvel,vector<vector<double> > &Press,vector<vector<double> > &zx,vector<vector<double> > &ey,vector<vector<double> > &ex,vector<vector<double> > &zy,double Re,double eps,double ep,vector<vector<double> > &rcs,vector<vector<double> > &rus,vector<vector<double> > &rvs, vector<vector<double> > &rho1, vector<vector<double> > &rho2)
{
   
    vector<vector<double> > U (N,vector<double>(N, 0));
    vector<vector<double> > V (N,vector<double>(N, 0));
    vector<vector<double> > g11 (N,vector<double>(N, 0));
    vector<vector<double> > g22 (N,vector<double>(N, 0));
    vector<vector<double> > g12 (N,vector<double>(N, 0));
    
    double JChx[N][N], JChy[N][N], zxhx[N][N], zxhy[N][N], zyhx[N][N], zyhy[N][N], exhx[N][N], exhy[N][N], eyhx[N][N], eyhy[N][N], g11hx[N][N], g12hx[N][N], g12hy[N][N],  g22hy[N][N], Estar1[N][N][3], Estar2[N][N][3], Estarv1[N][N][3], Estarv2[N][N][3], dEs1[N][N][3], dEs2[N][N][3], dEsv1[N][N][3], dEsv2[N][N][3];
   
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j  < N; ++j)
        {
            U[i][j]=0;
            V[i][j]=0;
            g11[i][j]=0;
            g12[i][j]=0;
            g22[i][j]=0;
            rcs[i][j]=0;
            rus[i][j]=0;
            rvs[i][j]=0;
            JChx[i][j]=0;
            JChy[i][j]=0;
            zxhx[i][j]=0;
            zyhx[i][j]=0;
            exhx[i][j]=0;
            eyhx[i][j]=0;
            zxhy[i][j]=0;
            zyhy[i][j]=0;
            exhy[i][j]=0;
            eyhy[i][j]=0;
            g11hx[i][j]=0;
            // g11hy[i][j]=0;
            g12hx[i][j]=0;
            g12hy[i][j]=0;
            // g22hx[i][j]=0;
            g22hy[i][j]=0;
            for (int k = 0; k < 3; ++k)
            {
                Estar1[i][j][k]=0;
                Estar2[i][j][k]=0;
                Estarv1[i][j][k]=0;
                Estarv2[i][j][k]=0;
                dEs1[i][j][k]=0;
                dEs2[i][j][k]=0;
                dEsv1[i][j][k]=0;
                dEsv2[i][j][k]=0;
            }
        }
    }
    
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            U[i][j]=(xvel[i][j]*zx[i][j]+yvel[i][j]*zy[i][j]);
            V[i][j]=(xvel[i][j]*ex[i][j]+yvel[i][j]*ey[i][j]);
        }
    }
    
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            Estar1[i][j][0]=U[i][j]/JC[i][j];
            Estar1[i][j][1]=(xvel[i][j]*U[i][j]+Press[i][j]*zx[i][j])/JC[i][j];
            Estar1[i][j][2]=(yvel[i][j]*U[i][j]+Press[i][j]*zy[i][j])/JC[i][j];
            Estar2[i][j][0]=V[i][j]/JC[i][j];
            Estar2[i][j][1]=(xvel[i][j]*V[i][j]+Press[i][j]*ex[i][j])/JC[i][j];
            Estar2[i][j][2]=(yvel[i][j]*V[i][j]+Press[i][j]*ey[i][j])/JC[i][j];
        }
    }
    
    for(int k=0;k<3; k++)
    {
        for (int i=1; i<N-1; i++)
        {
            for (int j=1; j<N-1;j++)
            {
                dEs1[i][j][k]=0.5*(Estar1[i][j+1][k]-Estar1[i][j-1][k]);
                dEs2[i][j][k]=0.5*(Estar2[i+1][j][k]-Estar2[i-1][j][k]);
            }
        }
    }
    
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            g11[i][j]=(pow(zx[i][j],2))+(pow(zy[i][j],2));
            g12[i][j]=zx[i][j]*ex[i][j]+zy[i][j]*ey[i][j];
            g22[i][j]=(pow(ex[i][j],2))+(pow(ey[i][j],2));
        }
    }
    for (int i=0; i<N-1; i++)
    {
        for (int j=0; j<N; j++)
        {
            JChy[i][j]=0.5*(JC[i+1][j]+JC[i][j]);
            zxhy[i][j]=0.5*(zx[i+1][j]+zx[i][j]);
            zyhy[i][j]=0.5*(zy[i+1][j]+zy[i][j]);
            exhy[i][j]=0.5*(ex[i+1][j]+ex[i][j]);
            eyhy[i][j]=0.5*(ey[i+1][j]+ey[i][j]);
        }
    }
    for (int j=0; j<N; j++)
    {
        JChy[N-1][j]=0.5*(JC[N-2][j]+JC[N-1][j]);
        zxhy[N-1][j]=0.5*(zx[N-2][j]+zx[N-1][j]);
        zyhy[N-1][j]=0.5*(zy[N-2][j]+zy[N-1][j]);
        exhy[N-1][j]=0.5*(ex[N-2][j]+ex[N-1][j]);
        eyhy[N-1][j]=0.5*(ey[N-2][j]+ey[N-1][j]);
    }
    for (int i=0; i<N; i++)
    { 
        for (int j=0; j<N-1; j++)
        {
            JChx[i][j]=0.5*(JC[i][j+1]+JC[i][j]);
            zxhx[i][j]=0.5*(zx[i][j+1]+zx[i][j]);
            zyhx[i][j]=0.5*(zy[i][j+1]+zy[i][j]);
            exhx[i][j]=0.5*(ex[i][j+1]+ex[i][j]);
            eyhx[i][j]=0.5*(ey[i][j+1]+ey[i][j]);
        }
    }
    for (int i=0; i<N; i++)
    {
        JChx[i][N-1]=0.5*(JC[i][N-2]+JC[i][N-1]);
        zyhx[i][N-1]=0.5*(zy[i][N-2]+zy[i][N-1]);
        zxhx[i][N-1]=0.5*(zx[i][N-2]+zx[i][N-1]);
        exhx[i][N-1]=0.5*(ex[i][N-2]+ex[i][N-1]);
        eyhx[i][N-1]=0.5*(ey[i][N-2]+ey[i][N-1]);
    }
    // #pragma omp parallel for schedule(guided,8)
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            g11hx[i][j]=(pow(zxhx[i][j],2))+(pow(zyhx[i][j],2));
            g12hx[i][j]=zxhx[i][j]*exhx[i][j]+zyhx[i][j]*eyhx[i][j];
            // g22hx[i][j]=(pow(exhx[i][j],2))+(pow(eyhx[i][j],2));
            // g11hy[i][j]=(pow(zxhy[i][j],2))+(pow(zyhy[i][j],2));
            g12hy[i][j]=zxhy[i][j]*exhy[i][j]+zyhy[i][j]*eyhy[i][j];
            g22hy[i][j]=(pow(exhy[i][j],2))+(pow(eyhy[i][j],2));
        }
    }
    for (int i=0; i<N-1; i++)
    {
        for (int j=0; j<N-1; j++)
        {
            Estarv1[i][j][0]=0;
             Estarv1[i][j][1]=(g11hx[i][j]*(xvel[i][j+1]-xvel[i][j])+g12hx[i][j]*(xvel[i+1][j]-xvel[i][j]))/(Re*JChx[i][j]);
            Estarv1[i][j][2]=(g11hx[i][j]*(yvel[i][j+1]-yvel[i][j])+g12hx[i][j]*(yvel[i+1][j]-yvel[i][j]))/(Re*JChx[i][j]);
        }
    }
    for (int i=0; i<N-1; i++)
    {
        for (int j=0; j<N-1; j++)
        {
            Estarv2[i][j][0]=0;
            Estarv2[i][j][1]=(g12hy[i][j]*(xvel[i][j+1]-xvel[i][j])+g22hy[i][j]*(xvel[i+1][j]-xvel[i][j]))/(Re*JChy[i][j]);
            Estarv2[i][j][2]=(g12hy[i][j]*(yvel[i][j+1]-yvel[i][j])+g22hy[i][j]*(yvel[i+1][j]-yvel[i][j]))/(Re*JChy[i][j]);
        }
    }
    for (int i=1; i<N-1; i++)
    {
        for (int j=1; j<N-1;j++)
        {
            dEsv1[i][j][0]=(Estarv1[i][j][0]-Estarv1[i][j-1][0]);
            dEsv1[i][j][1]=(Estarv1[i][j][1]-Estarv1[i][j-1][1]);
            dEsv1[i][j][2]=(Estarv1[i][j][2]-Estarv1[i][j-1][2]);
            dEsv2[i][j][0]=(Estarv2[i][j][0]-Estarv2[i-1][j][0]);
            dEsv2[i][j][1]=(Estarv2[i][j][1]-Estarv2[i-1][j][1]);
            dEsv2[i][j][2]=(Estarv2[i][j][2]-Estarv2[i-1][j][2]);
        }
    }
    // #pragma omp parallel for schedule(guided,8)
    for (int i=1; i<N-1; i++)
    {
        for (int j=1; j<N-1;j++)
        {
            rcs[i][j]=(dEs1[i][j][0]+dEs2[i][j][0]-dEsv1[i][j][0]-dEsv2[i][j][0])*JC[i][j];
            rus[i][j]=(dEs1[i][j][1]+dEs2[i][j][1]-dEsv1[i][j][1]-dEsv2[i][j][1])*JC[i][j];
            rvs[i][j]=(dEs1[i][j][2]+dEs2[i][j][2]-dEsv1[i][j][2]-dEsv2[i][j][2])*JC[i][j];
        }
    }
    spectralradius(N, JC, U, V, g11, g22,  rho1,  rho2);
    Diss( N, eps, rho1, rho2, Press, xvel, yvel, rcs, rus, rvs);
    smoothing( N, ep,  rcs,  rus,  rvs);
}
