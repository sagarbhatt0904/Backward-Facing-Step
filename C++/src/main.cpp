#include <iostream>
#include <cblas.h>
#include <cmath>
#include "GRIDGEN.hpp"
#include "INIT.hpp"
#include "METRIC.hpp"
#include "RHS.hpp"
#include <omp.h>
#include "contour.h"
#include <fstream>

using namespace std;

void BC (const int& N, vector<vector<double> > &u_new, vector<vector<double> > &v_new, vector<vector<double> > &p_new)
{
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < N; i++)
	{
		u_new[i][0] = 0.5 * (u_new[i][1] + u_new[i][N - 2]);
		u_new[0][i] = 0.5 * (u_new[1][i] + u_new[N - 2][i]);
		u_new[N - 1][i] = u_new[0][i];
		u_new[i][N - 1] = u_new[i][0];

		v_new[i][0] = 0.5 * (v_new[i][1] + v_new[i][N - 2]);
		v_new[0][i] = 0.5 * (v_new[1][i] + v_new[N - 2][i]);
		v_new[N - 1][i] = v_new[0][i];
		v_new[i][N - 1] = v_new[i][0];

		p_new[i][0] = 0.5 * (p_new[i][1] + p_new[i][N - 2]);
		p_new[0][i] = 0.5 * (p_new[1][i] + p_new[N - 2][i]);
		p_new[N - 1][i] = p_new[0][i];
		p_new[i][N - 1] = p_new[i][0];
	}

}

int main(int argc, char const *argv[])
{
	int N = 11; // Number of grid points
	int Re = 50;
	double grid_length = 10, grid_width = 2;
	std::vector<double> x;
	std::vector<double> y;

	vector<vector<double> > xvel (N, vector<double>(N, 0));
	vector<vector<double> > yvel (N, vector<double>(N, 0));
	vector<vector<double> > Press (N, vector<double>(N, 0));
	vector<vector<double> > J (N, vector<double>(N, 0));
	vector<vector<double> > ey (N, vector<double>(N, 0));
	vector<vector<double> > zx (N, vector<double>(N, 0));

	// Block to delete the objects since they are no longer necessary
	{
		GRIDGEN X(N, grid_length, x);
		GRIDGEN Y(N, grid_width, y);
		INIT init(N, x, y, xvel, yvel, Press);
		METRIC metric( N, x, y, zx, ey, J);
	}

	vector<vector<double> > u_new, u_new1, u_new2, u_new3;
	vector<vector<double> > v_new, v_new1, v_new2, v_new3;
	vector<vector<double> > p_new, p_new1, p_new2, p_new3;
	vector<vector<double> > U, V, g11, g22, rho1, rho2, dtau;
	Resize2D(N, rho1); Resize2D(N, rho2); Resize2D(N, U); Resize2D(N, V); Resize2D(N, g11); Resize2D(N, g22), Resize2D(N, dtau);
	Resize2D(N, u_new); Resize2D(N, u_new1); Resize2D(N, u_new2); Resize2D(N, u_new3);
	Resize2D(N, v_new); Resize2D(N, v_new1); Resize2D(N, v_new2); Resize2D(N, v_new3);
	Resize2D(N, p_new); Resize2D(N, p_new1); Resize2D(N, p_new2); Resize2D(N, p_new3);

	double norm = 1, CFL = 0.0001, VN = 0.01;
	int count, niter = 0;
	std::vector<double> normvec1D(N * N);
	char file1[] = "xvel.vtk", file2[] = "yvel.vtk", var1[] = "U", var2[] = "V";
	ofstream fout("norm.csv");
	while (norm >= pow(10, (-6)))
	{
		for (int i = 0; i < N; i++)
		{
			#pragma omp parallel for schedule(dynamic)
			for (int j = 0; j < N; j++)
			{
				U[i][j] = (xvel[i][j] * zx[i][j]);
				V[i][j] = (yvel[i][j] * ey[i][j]);
				g11[i][j] = (pow(zx[i][j], 2));
				g22[i][j] = (pow(ey[i][j], 2));
			}
		}
		for (int i = 0; i < N; i++)
		{
			#pragma omp parallel for schedule(dynamic)
			for (int j = 0; j < N; j++)
			{
				rho1[i][j] = (1 / J[i][j]) * (abs(U[i][j]) + sqrt((pow(U[i][j], 2) + g11[i][j])));
				rho2[i][j] = (1 / J[i][j]) * (abs(V[i][j]) + sqrt((pow(V[i][j], 2) + g22[i][j])));
			}
		}

		for (int i = 0; i < N; i++)
		{
			#pragma omp parallel for schedule(dynamic)
			for (int j = 0; j < N; j++)
			{
				rho1[i][j] = (1 / J[i][j]) * (abs(U[i][j]) + sqrt((pow(U[i][j], 2) + g11[i][j])));
				rho2[i][j] = (1 / J[i][j]) * (abs(V[i][j]) + sqrt((pow(V[i][j], 2) + g22[i][j])));
			}
		}

		for (int i = 0; i < N; ++i)
		{
			#pragma omp parallel for schedule(dynamic)
			for (int j = 0; j < N; ++j)
			{
				dtau[i][j] = min(CFL * (1 / J[i][j]) / (max(rho1[i][j], rho2[i][j])), (Re * VN / max(g11[i][j], g22[i][j])));
			}
		}
		/*Fourth order Runge-Kutta*/

		RHS step1(N, Re, xvel, yvel, Press, zx, ey, J);

		for (int i = 1; i < N - 1; i++)			// First step of RK
		{
			#pragma omp parallel for schedule(dynamic)
			for (int j = 1; j < N - 1; j++)
			{
				p_new1[i][j] = Press[i][j] + 0.25 * (dtau[i][j] * step1.rhs[0][i][j]);
				u_new1[i][j] = xvel[i][j] + 0.25 * (dtau[i][j] * step1.rhs[1][i][j]);
				v_new1[i][j] = yvel[i][j] + 0.25 * (dtau[i][j] * step1.rhs[2][i][j]);
			}
		}
		BC(N, u_new1, v_new1, p_new1);			// BC after first step RK

		RHS step2(N, Re, u_new1, v_new1, p_new1, zx, ey, J);

		for (int i = 1; i < N - 1; i++)			// Second step of RK
		{
			#pragma omp parallel for schedule(dynamic)
			for (int j = 1; j < N - 1; j++)
			{
				p_new2[i][j] = Press[i][j] + 0.33 * (dtau[i][j] * step2.rhs[0][i][j]);
				u_new2[i][j] = xvel[i][j] + 0.33 * (dtau[i][j] * step2.rhs[1][i][j]);
				v_new2[i][j] = yvel[i][j] + 0.33 * (dtau[i][j] * step2.rhs[2][i][j]);
			}
		}
		BC(N, u_new2, v_new2, p_new2);			// BC after second step RK

		RHS step3(N, Re, u_new2, v_new2, p_new2, zx, ey, J);

		for (int i = 1; i < N - 1; i++)			// Third step of RK
		{
			#pragma omp parallel for schedule(dynamic)
			for (int j = 1; j < N - 1; j++)
			{

				p_new3[i][j] = Press[i][j] + 0.5 * (dtau[i][j] * step3.rhs[0][i][j]);
				u_new3[i][j] = xvel[i][j] + 0.5 * (dtau[i][j] * step3.rhs[1][i][j]);
				v_new3[i][j] = yvel[i][j] + 0.5 * (dtau[i][j] * step3.rhs[2][i][j]);
			}
		}
		BC(N, u_new3, v_new3, p_new3);			// BC after third step RK

		RHS step4(N, Re, u_new3, v_new3, p_new3, zx, ey, J);

		for (int i = 1; i < N - 1; i++)			// Fourth step of RK
		{
			#pragma omp parallel for schedule(dynamic)
			for (int j = 1; j < N - 1; j++)
			{
				p_new[i][j] = Press[i][j] + (dtau[i][j] * step4.rhs[0][i][j]);
				u_new[i][j] = xvel[i][j] + (dtau[i][j] * step4.rhs[1][i][j]);
				v_new[i][j] = yvel[i][j] + (dtau[i][j] * step4.rhs[2][i][j]);
			}
		}
		BC(N, u_new, v_new, p_new);			// BC after fourth step RK

		count = 0;
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				normvec1D[count] = abs(u_new[i][j] - xvel[i][j]);
				count++;
			}
		}

		norm = cblas_dnrm2(N * N, normvec1D.data(), 1);

		xvel = u_new;
		yvel = v_new;
		Press = p_new;
		niter++;
		if ((niter % 1000) == 0)
		{
			contour(u_new, N, file1, var1);
			contour(v_new, N, file2, var2);
			cout << norm << " " << niter << endl;
			fout << niter << "," << norm << endl;
		}
		// cout<<dtau[2][2]<<" "<<dtau[20][20]<<endl;

	}

	fout.close();
	return 0;
}