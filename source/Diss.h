#include<math.h>


// computing dissipation

void Diss(int N,double eps,double** rho1,double** rho2,double** Press,double** xvel,double** yvel,double** rcs,double** rus,double** rvs)

{ double D1_c[N][N], D2_c[N][N],D1_u[N][N], D2_u[N][N],D1_v[N][N], D2_v[N][N], Diss_c[N][N],Diss_u[N][N],Diss_v[N][N];

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j  < N; ++j)
        {
            D1_c[i][j]=0;
            D1_u[i][j]=0;
            D1_v[i][j]=0;
            D2_c[i][j]=0;
            D2_u[i][j]=0;
            D2_v[i][j]=0;
            Diss_c[i][j]=0;
            Diss_u[i][j]=0;
            Diss_v[i][j]=0;
        }
    }


    for (int i=2; i<N-2; i++)
    	{
    	    for (int j=2; j<N-2; j++)
        	    {
    		D1_c[i][j]=eps*(rho1[i][j]+rho1[i+1][j])*0.5*(Press[i+2][j]-3*Press[i+1][j]+3*Press[i][j]-Press[i-1][j]);
    		D1_u[i][j]=eps*(rho1[i][j]+rho1[i+1][j])*0.5*(xvel[i+2][j]-3*xvel[i+1][j]+3*xvel[i][j]-xvel[i-1][j]);
    		D1_v[i][j]=eps*(rho1[i][j]+rho1[i+1][j])*0.5*(yvel[i+2][j]-3*yvel[i+1][j]+3*yvel[i][j]-yvel[i-1][j]);
    		D2_c[i][j]=eps*(rho2[i][j]+rho2[i][j+1])*0.5*(Press[i][j+2]-3*Press[i][j+1]+3*Press[i][j]-Press[i][j-1]);
    		D2_u[i][j]=eps*(rho2[i][j]+rho2[i][j+1])*0.5*(xvel[i][j+2]-3*xvel[i][j+1]+3*xvel[i][j]-xvel[i][j-1]);
    		D2_v[i][j]=eps*(rho2[i][j]+rho2[i][j+1])*0.5*(yvel[i][j+2]-3*yvel[i][j+1]+3*yvel[i][j]-yvel[i][j-1]);
    		Diss_c[i][j]=((D1_c[i+1][j]+D1_c[i][j])*0.5-0.5*(D1_c[i][j]+D1_c[i-1][j]))+((D2_c[i+1][j]+D2_c[i][j])*0.5-0.5*(D2_c[i][j]+D2_c[i][j-1]));
    		Diss_u[i][j]=((D1_u[i+1][j]+D1_u[i][j])*0.5-0.5*(D1_u[i][j]+D1_u[i-1][j]))+((D2_u[i+1][j]+D2_u[i][j])*0.5-0.5*(D2_u[i][j]+D2_u[i][j-1]));
    		Diss_v[i][j]=((D1_v[i+1][j]+D1_v[i][j])*0.5-0.5*(D1_v[i][j]+D1_v[i-1][j]))+((D2_v[i+1][j]+D2_v[i][j])*0.5-0.5*(D2_v[i][j]+D2_v[i][j-1]));
            
        	     }
    	}
    for (int i=0; i<N; i++)
    {
        Diss_c[i][0]=Diss_c[i][1];
        Diss_u[i][0]=Diss_u[i][1];
        Diss_v[i][0]=Diss_v[i][1];
        Diss_c[i][N-2]=Diss_c[i][N-3];
        Diss_u[i][N-2]=Diss_u[i][N-3];
        Diss_v[i][N-2]=Diss_v[i][N-3];
        Diss_c[i][N-1]=Diss_c[i][N-2];
        Diss_u[i][N-1]=Diss_u[i][N-2];
        Diss_v[i][N-1]=Diss_v[i][N-2];
    }
    for (int j=0; j<N; j++)
    {
        Diss_c[0][j]=(Diss_c[1][j]);
        Diss_u[0][j]=(Diss_u[1][j]);
        Diss_v[0][j]=(Diss_v[1][j]);
        Diss_c[N-2][j]=Diss_c[N-3][j];
        Diss_u[N-2][j]=Diss_u[N-3][j];
        Diss_v[N-2][j]=Diss_v[N-3][j];
        Diss_c[N-1][j]=Diss_c[N-2][j];
        Diss_u[N-1][j]=Diss_u[N-2][j];
        Diss_v[N-1][j]=Diss_v[N-2][j];
    }
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            rcs[i][j]=rcs[i][j]+Diss_c[i][j];
            rus[i][j]=rus[i][j]+Diss_u[i][j];
            rvs[i][j]=rvs[i][j]+Diss_v[i][j];
        }
    }
}
