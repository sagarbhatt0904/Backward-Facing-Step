/****Tyler Green Vortex****
***************************

***************************/
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<omp.h>
#include"init.h"
#include"gridgen.h"
#include"metric.h"
#include"RHS.h"
#include"BC.h"
#include"timer.h"

using namespace std;

int main()
{
	timespec before, after, time_diff;
	int N= 61;                          
	double st=1; 			  
	double Re=50;
	double eps=0.01;
	double ep=0.01;
	double xmax=2;                         
	double ymax=20;

	get_time(&before);	
	double** x; double**y; double** xvel; double** yvel; double** xvel1; double** yvel1; double** Press; double** dummyu;
	x=new double* [N];  y=new double* [N];  xvel=new double* [N];  yvel=new double* [N];  xvel1=new double* [N];  yvel1=new double* [N];  Press=new double* [N]; dummyu=new double* [N];		  	
	
	#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
		x[i]=new double [N];  y[i]=new double [N];  xvel[i]=new double [N];  yvel[i]=new double [N];  xvel1[i]=new double [N];  yvel1[i]=new double [N];  Press[i]=new double [N]; dummyu[i]=new double [N];		  	
	}
	double norms, sum, norm_diff[N][N];    
	gridgen(N, st, xmax, ymax, x, y);  // Grid generation
	
	
	init( N, x, y, xvel, xvel1, yvel, yvel1, Press,dummyu); 
	
	// Initial conditions
	
	
	double** u_new; double** u_new1; double** u_new2; double** u_new3; double** v_new; double** v_new1; double** v_new2; double** v_new3; double** p_new; double** p_new1; double** p_new2; double** p_new3;
	double** u_k; double** u_k1; double** u_k2; double** u_k3; double** v_k; double** v_k1; double** v_k2; double** v_k3; double** p_k; double** p_k1; double** p_k2; double** p_k3;
	double** u_old; double** u_old1; double** u_old2; double** u_old3; double** v_old; double** v_old1; double** v_old2; double** v_old3; double** u_n; double** u_n1; double** u_n2; double** u_n3;
	double** v_n; double** v_n1; double** v_n2; double** v_n3; double** p_n; double** dtau; double** JC; double** ex; double** ey; double** zx; double** zy; double** rus; double** rvs; double** rcs;
	double** p1; double** p2; double** p3; double** rho1; double** rho2;

	u_new=new double* [N];  u_new1=new double* [N];  u_new2=new double* [N];  u_new3=new double* [N];  v_new=new double* [N];  v_new1=new double* [N];  v_new2=new double* [N];  v_new3=new double* [N];  p_new=new double* [N];  p_new1=new double* [N];  p_new2=new double* [N];  p_new3=new double* [N];
	u_k=new double* [N];  u_k1=new double* [N];  u_k2=new double* [N];  u_k3=new double* [N];  v_k=new double* [N];  v_k1=new double* [N];  v_k2=new double* [N];  v_k3=new double* [N];  p_k=new double* [N];  p_k1=new double* [N];  p_k2=new double* [N];  p_k3=new double* [N];
	u_old=new double* [N];  u_old1=new double* [N];  u_old2=new double* [N];  u_old3=new double* [N];  v_old=new double* [N];  v_old1=new double* [N];  v_old2=new double* [N];  v_old3=new double* [N];  u_n=new double* [N];  u_n1=new double* [N];  u_n2=new double* [N];  u_n3=new double* [N];
	v_n=new double* [N];  v_n1=new double* [N];  v_n2=new double* [N];  v_n3=new double* [N];  p_n=new double* [N];  dtau=new double* [N];  JC=new double* [N];  ex=new double* [N];  ey=new double* [N];  zx=new double* [N];  zy=new double* [N];  rus=new double* [N];  rvs=new double* [N];  rcs=new double* [N];
	p1=new double* [N];  p2=new double* [N];  p3=new double* [N];rho1=new double* [N];rho2=new double* [N];

	#pragma omp parallel for
	for (int i = 0; i < N; ++i)
	{
		 u_new[i]=new double [N];  u_new1[i]=new double [N];  u_new2[i]=new double [N];  u_new3[i]=new double [N];  v_new[i]=new double [N];  v_new1[i]=new double [N];  v_new2[i]=new double [N];  v_new3[i]=new double [N];  p_new[i]=new double [N];  p_new1[i]=new double [N];  p_new2[i]=new double [N];  p_new3[i]=new double [N];
		 u_k[i]=new double [N];  u_k1[i]=new double [N];  u_k2[i]=new double [N];  u_k3[i]=new double [N];  v_k[i]=new double [N];  v_k1[i]=new double [N];  v_k2[i]=new double [N];  v_k3[i]=new double [N];  p_k[i]=new double [N];  p_k1[i]=new double [N];  p_k2[i]=new double [N];  p_k3[i]=new double [N];
		 u_old[i]=new double [N];  u_old1[i]=new double [N];  u_old2[i]=new double [N];  u_old3[i]=new double [N];  v_old[i]=new double [N];  v_old1[i]=new double [N];  v_old2[i]=new double [N];  v_old3[i]=new double [N];  u_n[i]=new double [N];  u_n1[i]=new double [N];  u_n2[i]=new double [N];  u_n3[i]=new double [N];
		 v_n[i]=new double [N];  v_n1[i]=new double [N];  v_n2[i]=new double [N];  v_n3[i]=new double [N];  p_n[i]=new double [N];  dtau[i]=new double [N];  JC[i]=new double [N];  ex[i]=new double [N];  ey[i]=new double [N];  zx[i]=new double [N];  zy[i]=new double [N];  rus[i]=new double [N];  rvs[i]=new double [N];  rcs[i]=new double [N];
		  p1[i]=new double [N];  p2[i]=new double [N];  p3[i]=new double [N]; rho1[i]=new double [N];rho2[i]=new double [N];
	}

	metric( N, x, y,zx,zy,ex,ey,JC);   //Metric Calculation
	
	#pragma omp parallel for
	for (int i=0; i<N; i++)  // Initializing all the velocity variables and dtau
	{
		for (int j=0; j<N; j++)
		{
		
			u_new1[i][j]=xvel[i][j];
			u_new2[i][j]=xvel[i][j];
			u_new3[i][j]=xvel[i][j];
			u_new[i][j]=xvel[i][j];
			u_k1[i][j]=xvel[i][j];
			u_k2[i][j]=xvel[i][j];
			u_k3[i][j]=xvel[i][j];
			u_k[i][j]=xvel[i][j];
			u_n1[i][j]=xvel[i][j];
			u_n2[i][j]=xvel[i][j];
			u_n3[i][j]=xvel[i][j];
			u_n[i][j]=xvel[i][j];
			u_old1[i][j]=xvel[i][j];
			u_old2[i][j]=xvel[i][j];
			u_old3[i][j]=xvel[i][j];
			u_old[i][j]=xvel[i][j];
			v_new1[i][j]=yvel[i][j];
			v_new2[i][j]=yvel[i][j];
			v_new3[i][j]=yvel[i][j];
			v_new[i][j]=yvel[i][j];
			v_k1[i][j]=yvel[i][j];
			v_k2[i][j]=yvel[i][j];
			v_k3[i][j]=yvel[i][j];
			v_k[i][j]=yvel[i][j];
			v_n1[i][j]=yvel[i][j];
			v_n2[i][j]=yvel[i][j];
			v_n3[i][j]=yvel[i][j];
			v_n[i][j]=yvel[i][j];
			v_old1[i][j]=yvel[i][j];
			v_old2[i][j]=yvel[i][j];
			v_old3[i][j]=yvel[i][j];
			v_old[i][j]=yvel[i][j];
			p_new1[i][j]=Press[i][j];
			p_new2[i][j]=Press[i][j];
			p_new3[i][j]=Press[i][j];
			p_new[i][j]=Press[i][j];
			p_n[i][j]=Press[i][j];
		    dtau[i][j]=0.0001;
		    	
		}
	}
	// double dj=1;

	//while( dj>=0.00001)
	for (int lop = 0; lop < 100; ++lop)
	{
	 	// sum=0;   

	 

	 	/*Fourth order Runge-Kutta*/

	    RHS(N,JC,u_k,v_k,Press,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
	    
	    
		#pragma omp parallel for
		for (int i =1; i<N-1; i++)			// First step of RK
	    {    for (int j =1; j<N-1; j++)
	        {
	               	    
	            p_new1[i][j]=Press[i][j]+0.25*(dtau[i][j]*rcs[i][j]);
	            u_new1[i][j]=u_k[i][j]+0.25*(dtau[i][j]*rus[i][j]);
	            v_new1[i][j]=v_k[i][j]+0.25*(dtau[i][j]*rvs[i][j]);
	        }
	    }
	    BC(N,u_new1,v_new1,p_new1,dummyu);			// BC after first step RK
	    
	    RHS(N,JC,u_new1,v_new1,p_new1,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
	    
		#pragma omp parallel for
	    for (int i =1; i<N-1; i++)			// Second step of RK
	    {
	        for (int j =1; j<N-1; j++)
	        {
	            p_new2[i][j]=Press[i][j]+0.33*(dtau[i][j]*rcs[i][j]);
	            u_new2[i][j]=u_k[i][j]+0.33*(dtau[i][j]*rus[i][j]);
	            v_new2[i][j]=v_k[i][j]+0.33*(dtau[i][j]*rvs[i][j]);
	        }
	    }
	    BC(N,u_new2,v_new2,p_new2,dummyu);			// BC after second step RK
	    
	    RHS(N,JC,u_new2,v_new2,p_new2,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
	    
		#pragma omp parallel for
	    for (int i =1; i<N-1; i++)			// Third step of RK
	    {
	        for (int j =1; j<N-1; j++)
	        {   
	            
	            p_new3[i][j]=Press[i][j]+0.5*(dtau[i][j]*rcs[i][j]);
	            u_new3[i][j]=u_k[i][j]+0.5*(dtau[i][j]*rus[i][j]);
	            v_new3[i][j]=v_k[i][j]+0.5*(dtau[i][j]*rvs[i][j]);
	        }
	    }
	    BC(N,u_new3,v_new3,p_new3,dummyu);			// BC after third step RK
	    
	    RHS(N,JC,u_new3,v_new3,p_new3,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
	    
		#pragma omp parallel for
	    for (int i =1; i<N-1; i++)			// Fourth step of RK
	    {
	        for (int j =1; j<N-1; j++)
	        {   
	           
	            p_new[i][j]=Press[i][j]+(dtau[i][j]*rcs[i][j]);
	            u_new[i][j]=u_k[i][j]+(dtau[i][j]*rus[i][j]);
	            v_new[i][j]=v_k[i][j]+(dtau[i][j]*rvs[i][j]);
	        }
	    }
	    BC(N,u_new,v_new,p_new,dummyu);			// BC after fourth step RK`

		#pragma omp parallel for
	    for (int i =0; i<N; i++)			// Updating old values
	    {
	        for (int j =0; j<N; j++)
	        { 
			    u_old1[i][j]=u_n1[i][j];
			    u_n1[i][j]=u_k1[i][j];
			    u_k1[i][j]=u_new1[i][j];
			    u_old2[i][j]=u_n2[i][j];
			    u_n2[i][j]=u_k2[i][j];
			    u_k2[i][j]=u_new2[i][j];
			    u_old3[i][j]=u_n3[i][j];
			    u_n3[i][j]=u_k3[i][j];
			    u_k3[i][j]=u_new3[i][j];
			    u_old[i][j]=u_n[i][j];
			    u_n[i][j]=u_k[i][j];
			    u_k[i][j]=u_new[i][j];
			    v_old1[i][j]=v_n1[i][j];
			    v_n1[i][j]=v_k1[i][j];
			    v_k1[i][j]=v_new1[i][j];
			    v_old2[i][j]=v_n2[i][j];
			    v_n2[i][j]=v_k2[i][j];
			    v_k2[i][j]=v_new2[i][j];
			    v_old3[i][j]=v_n3[i][j];
			    v_n3[i][j]=v_k3[i][j];
			    v_k3[i][j]=v_new3[i][j];
			    v_old[i][j]=v_n[i][j];
			    v_n[i][j]=v_k[i][j];
			    v_k[i][j]=v_new[i][j];		    
			    p1[i][j]=p_new1[i][j];
			    p2[i][j]=p_new2[i][j];
			    p3[i][j]=p_new3[i][j];
			    Press[i][j]=p_new[i][j];
		  }
	    } 		   
	    
	    
	}  	
	get_time(&after);
	diff(&before,&after,&time_diff);
	
	double time_s = time_diff.tv_sec + (double)(time_diff.tv_nsec)/1.0e9;
	cout<<"Time: "<<time_s<<endl;
    /* Writing Data to file */
	    ofstream xout("xvel.vtk");
		xout<<"# vtk DataFile Version 2.0\nVTK from matlab\n"<<"ASCII\n"<<"DATASET STRUCTURED_POINTS\n"<<"DIMENSIONS "<<N<<" "<<N<<" "<<1<<"\nSPACING 1 1 1 \n"<<"ORIGIN 0 0 0\n"<<"POINT_DATA "<<N*N<<"\nSCALARS xvel float 1\n"<<"LOOKUP_TABLE default\n";
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				xout<<u_new[i][j]<<" ";
			}
		}
		
		xout.close();
		
		ofstream yout("yvel.vtk");
		yout<<"# vtk DataFile Version 2.0\nVTK from matlab\n"<<"ASCII\n"<<"DATASET STRUCTURED_POINTS\n"<<"DIMENSIONS "<<N<<" "<<N<<" "<<1<<"\nSPACING 1 1 1 \n"<<"ORIGIN 0 0 0\n"<<"POINT_DATA "<<N*N<<"\nSCALARS yvel float 1\n"<<"LOOKUP_TABLE default\n";
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				yout<<v_new[i][j]<<" ";
			}
		}
		yout.close();
	
	
    return 0;
}
