/****Tyler Green Vortex****
***************************

***************************/
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <vector>
#include "init.h"
#include "gridgen.h"
#include "metric.h"
#include "RHS.h"
#include "BC.h"

using namespace std;

int main()
{
	
	int N= 61;                          
	double st=1; 			  
	double Re=50;
	double eps=0.01;
	double ep=0.01;
	double xmax=20;                         
	double ymax=2;

	// get_time(&before);
	vector<vector<double> > x (N,vector<double>(N, 0));
	vector<vector<double> > y (N,vector<double>(N, 0));
	vector<vector<double> > xvel (N,vector<double>(N, 0));	
	vector<vector<double> > yvel (N,vector<double>(N, 0));
	vector<vector<double> > xvel1 (N,vector<double>(N, 0));
	vector<vector<double> > yvel1 (N,vector<double>(N, 0));
	vector<vector<double> > Press (N,vector<double>(N, 0));
	vector<vector<double> > dummyu (N,vector<double>(N, 0));
	
	
	gridgen(N, st, xmax, ymax, x, y);  // Grid generation
	
	
	init( N, x, y, xvel, xvel1, yvel, yvel1, Press,dummyu); // Initial conditions
	
	vector<vector<double> > u_new (N,vector<double>(N, 0));
	vector<vector<double> > u_new1 (N,vector<double>(N, 0));
	vector<vector<double> > u_new2 (N,vector<double>(N, 0));
	vector<vector<double> > u_new3 (N,vector<double>(N, 0));
	vector<vector<double> > u_k (N,vector<double>(N, 0));
	vector<vector<double> > u_k1 (N,vector<double>(N, 0));
	vector<vector<double> > u_k2 (N,vector<double>(N, 0));
	vector<vector<double> > u_k3 (N,vector<double>(N, 0));
	vector<vector<double> > u_n (N,vector<double>(N, 0));
	vector<vector<double> > u_n1 (N,vector<double>(N, 0));
	vector<vector<double> > u_n2 (N,vector<double>(N, 0));
	vector<vector<double> > u_n3 (N,vector<double>(N, 0));
	vector<vector<double> > u_old (N,vector<double>(N, 0));
	vector<vector<double> > u_old1 (N,vector<double>(N, 0));
	vector<vector<double> > u_old2 (N,vector<double>(N, 0));
	vector<vector<double> > u_old3 (N,vector<double>(N, 0));

	vector<vector<double> > v_new (N,vector<double>(N, 0));
	vector<vector<double> > v_new1 (N,vector<double>(N, 0));
	vector<vector<double> > v_new2 (N,vector<double>(N, 0));
	vector<vector<double> > v_new3 (N,vector<double>(N, 0));
	vector<vector<double> > v_k (N,vector<double>(N, 0));
	vector<vector<double> > v_k1 (N,vector<double>(N, 0));
	vector<vector<double> > v_k2 (N,vector<double>(N, 0));
	vector<vector<double> > v_k3 (N,vector<double>(N, 0));
	vector<vector<double> > v_n (N,vector<double>(N, 0));
	vector<vector<double> > v_n1 (N,vector<double>(N, 0));
	vector<vector<double> > v_n2 (N,vector<double>(N, 0));
	vector<vector<double> > v_n3 (N,vector<double>(N, 0));
	vector<vector<double> > v_old (N,vector<double>(N, 0));
	vector<vector<double> > v_old1 (N,vector<double>(N, 0));
	vector<vector<double> > v_old2 (N,vector<double>(N, 0));
	vector<vector<double> > v_old3 (N,vector<double>(N, 0));

	vector<vector<double> > p_new (N,vector<double>(N, 0));
	vector<vector<double> > p_new1 (N,vector<double>(N, 0));
	vector<vector<double> > p_new2 (N,vector<double>(N, 0));
	vector<vector<double> > p_new3 (N,vector<double>(N, 0));
	vector<vector<double> > p_k (N,vector<double>(N, 0));
	vector<vector<double> > p_k1 (N,vector<double>(N, 0));
	vector<vector<double> > p_k2 (N,vector<double>(N, 0));
	vector<vector<double> > p_k3 (N,vector<double>(N, 0));
	vector<vector<double> > p_n (N,vector<double>(N, 0));
	vector<vector<double> > p_n1 (N,vector<double>(N, 0));
	vector<vector<double> > p_n2 (N,vector<double>(N, 0));
	vector<vector<double> > p_n3 (N,vector<double>(N, 0));
	vector<vector<double> > p_old (N,vector<double>(N, 0));
	vector<vector<double> > p_old1 (N,vector<double>(N, 0));
	vector<vector<double> > p_old2 (N,vector<double>(N, 0));
	vector<vector<double> > p_old3 (N,vector<double>(N, 0));

	vector<vector<double> > dtau (N,vector<double>(N, 0));
	vector<vector<double> > JC (N,vector<double>(N, 0));
	vector<vector<double> > ex (N,vector<double>(N, 0));
	vector<vector<double> > ey (N,vector<double>(N, 0));
	vector<vector<double> > zx (N,vector<double>(N, 0));
	vector<vector<double> > zy (N,vector<double>(N, 0));
	vector<vector<double> > rus (N,vector<double>(N, 0));
	vector<vector<double> > rvs (N,vector<double>(N, 0));
	vector<vector<double> > rcs (N,vector<double>(N, 0));
	vector<vector<double> > rho1 (N,vector<double>(N, 0));
	vector<vector<double> > rho2(N,vector<double>(N, 0));


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
	for (int lop = 0; lop < 100000; ++lop)
	{
	   

	 

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
			    Press[i][j]=p_new[i][j];
		  }
	    } 		   
	    
	    
	}  	
	
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
