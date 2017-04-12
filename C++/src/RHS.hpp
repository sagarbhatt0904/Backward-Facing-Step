#include <vector>
#include <cmath>
#include <omp.h>
#include "TRI.h"

/*

Note: Improve this!!
Spatially atleast!

*/

void Resize2D(const int& N, std::vector<std::vector<double> >& v)
{
    v.resize(N);
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; ++i)
    {
        v[i].resize(N);
    }
}

class CONVECTIVE
{

public:
    std::vector<std::vector<double> > U;
    std::vector<std::vector<double> > V;
    std::vector<std::vector<std::vector<double> > > Estar1;
    std::vector<std::vector<std::vector<double> > > Estar2;
    CONVECTIVE() {}
    CONVECTIVE(const int &,
               const std::vector<std::vector<double> >& ,
               const std::vector<std::vector<double> >& ,
               const std::vector<std::vector<double> >& ,
               const std::vector<std::vector<double> >& ,
               const std::vector<std::vector<double> >& ,
               const std::vector<std::vector<double> >& );


};

CONVECTIVE::CONVECTIVE(const int &N,
                       const std::vector<std::vector<double> >& xvel,
                       const std::vector<std::vector<double> >& yvel,
                       const std::vector<std::vector<double> >& Press,
                       const std::vector<std::vector<double> >& zx,
                       const std::vector<std::vector<double> >& ey,
                       const std::vector<std::vector<double> >& J)
{

    Resize2D(N, U); Resize2D(N, V);
    Estar1.resize(3);
    Estar2.resize(3);

    for (int i = 0; i < 3; ++i)
    {
        Estar1[i].resize(N);
        Estar2[i].resize(N);
    }
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Estar1[j][i].resize(N);
            Estar2[j][i].resize(N);
        }
    }
    for (int i = 0; i < N; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++)
        {
            U[i][j] = (xvel[i][j] * zx[i][j]);
            V[i][j] = (yvel[i][j] * ey[i][j]);
        }
    }

    for (int i = 0; i < N; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++)
        {
            Estar1[0][i][j] = U[i][j] / J[i][j];
            Estar1[1][i][j] = (xvel[i][j] * U[i][j] + Press[i][j] * zx[i][j]) / J[i][j];
            Estar1[2][i][j] = (yvel[i][j] * U[i][j]) / J[i][j];
            Estar2[0][i][j] = V[i][j] / J[i][j];
            Estar2[1][i][j] = (xvel[i][j] * V[i][j]) / J[i][j];
            Estar2[2][i][j] = (yvel[i][j] * V[i][j] + Press[i][j] * ey[i][j]) / J[i][j];
        }
    }
}


class VISCOUS
{
    std::vector<std::vector<double> > Jhy;
    std::vector<std::vector<double> > Jhx;
    std::vector<std::vector<double> > zxhx;
    std::vector<std::vector<double> > eyhy;
    std::vector<std::vector<double> > zxhy;
    std::vector<std::vector<double> > eyhx;
    std::vector<std::vector<double> > g11hx;
    std::vector<std::vector<double> > g22hy;

public:
    std::vector<std::vector<double> > g11;
    std::vector<std::vector<double> > g22;
    std::vector<std::vector<std::vector<double> > > Estarv1;
    std::vector<std::vector<std::vector<double> > > Estarv2;
    VISCOUS() {}
    VISCOUS(   const int &,
               const int &,
               const std::vector<std::vector<double> >& ,
               const std::vector<std::vector<double> >& ,
               const std::vector<std::vector<double> >& ,
               const std::vector<std::vector<double> >& ,
               const std::vector<std::vector<double> >& );

};

VISCOUS::VISCOUS  (const int &N,
                   const int &Re,
                   const std::vector<std::vector<double> >& xvel,
                   const std::vector<std::vector<double> >& yvel,
                   const std::vector<std::vector<double> >& zx,
                   const std::vector<std::vector<double> >& ey,
                   const std::vector<std::vector<double> >& J)
{

    Estarv1.resize(3);
    Estarv2.resize(3);
    Resize2D(N, Jhy); Resize2D(N, Jhx); Resize2D(N, zxhx); Resize2D(N, zxhy); Resize2D(N, eyhx); Resize2D(N, eyhy); Resize2D(N, g11); Resize2D(N, g22);
    Resize2D(N, g11hx); Resize2D(N, g22hy);

    for (int i = 0; i < 3; ++i)
    {
        Estarv1[i].resize(N);
        Estarv2[i].resize(N);
    }
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Estarv1[j][i].resize(N);
            Estarv2[j][i].resize(N);
        }
    }
    for (int i = 0; i < N - 1; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++)
        {
            Jhy[i][j] = 0.5 * (J[i + 1][j] + J[i][j]);
            zxhy[i][j] = 0.5 * (zx[i + 1][j] + zx[i][j]);
            eyhy[i][j] = 0.5 * (ey[i + 1][j] + ey[i][j]);

        }
    }
    for (int i = 0; i < N; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N - 1; j++)
        {
            Jhx[i][j] = 0.5 * (J[i][j + 1] + J[i][j]);
            zxhx[i][j] = 0.5 * (zx[i][j + 1] + zx[i][j]);
            eyhx[i][j] = 0.5 * (ey[i][j + 1] + ey[i][j]);
        }
    }

    #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < N; j++)
    {
        Jhy[N - 1][j] = 0.5 * (J[N - 2][j] + J[N - 1][j]);
        zxhy[N - 1][j] = 0.5 * (zx[N - 2][j] + zx[N - 1][j]);
        eyhy[N - 1][j] = 0.5 * (ey[N - 2][j] + ey[N - 1][j]);
        Jhx[j][N - 1] = 0.5 * (J[j][N - 2] + J[j][N - 1]);
        zxhx[j][N - 1] = 0.5 * (zx[j][N - 2] + zx[j][N - 1]);
        eyhx[j][N - 1] = 0.5 * (ey[j][N - 2] + ey[j][N - 1]);
    }


    for (int i = 0; i < N; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++)
        {
            g11[i][j] = (pow(zx[i][j], 2));
            g22[i][j] = (pow(ey[i][j], 2));
            g11hx[i][j] = (pow(zxhx[i][j], 2));
            g22hy[i][j] = (pow(eyhy[i][j], 2));
        }
    }
    for (int i = 0; i < N; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N - 1; j++)
        {
            Estarv1[0][i][j] = 0;
            Estarv1[1][i][j] = (g11hx[i][j] * (xvel[i][j + 1] - xvel[i][j])) / (Re * Jhx[i][j]);
            Estarv1[2][i][j] = (g11hx[i][j] * (yvel[i][j + 1] - yvel[i][j])) / (Re * Jhx[i][j]);
        }
    }
    for (int i = 0; i < N - 1; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++)
        {
            Estarv2[0][i][j] = 0;
            Estarv2[1][i][j] = (g22hy[i][j] * (xvel[i + 1][j] - xvel[i][j])) / (Re * Jhy[i][j]);
            Estarv2[2][i][j] = (g22hy[i][j] * (yvel[i + 1][j] - yvel[i][j])) / (Re * Jhy[i][j]);
        }

    }
}


class RHS: public CONVECTIVE, public VISCOUS
{

    std::vector<std::vector<double> > D1_c, D2_c, D1_u, D2_u, D1_v, D2_v, Diss_c, Diss_u , Diss_v, rho1, rho2   ;
    double eps, ep;
public:
    std::vector<std::vector<std::vector<double> > > rhs;
    RHS();
    RHS(    const int &,
            const int &,
            const std::vector<std::vector<double> >& ,
            const std::vector<std::vector<double> >& ,
            const std::vector<std::vector<double> >& ,
            const std::vector<std::vector<double> >& ,
            const std::vector<std::vector<double> >& ,
            const std::vector<std::vector<double> >& );
    void spectralradius(    const int& ,
                            const vector<vector<double> > & ,
                            const vector<vector<double> > & ,
                            const vector<vector<double> > & ,
                            const vector<vector<double> > & ,
                            const vector<vector<double> > & );
    void Diss(  const int& ,
                const vector<vector<double> > &,
                const vector<vector<double> > &,
                const vector<vector<double> > &);
    void smoothing(const int&);

};

RHS::RHS(   const int &N,
            const int &Re,
            const std::vector<std::vector<double> >& xvel,
            const std::vector<std::vector<double> >& yvel,
            const std::vector<std::vector<double> >& Press,
            const std::vector<std::vector<double> >& zx,
            const std::vector<std::vector<double> >& ey,
            const std::vector<std::vector<double> >& J)
{

    rhs.resize(3);

    for (int i = 0; i < 3; ++i)
    {
        rhs[i].resize(N);
    }
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            rhs[j][i].resize(N);
        }
    }
    CONVECTIVE convective(N, xvel, yvel, Press, zx, ey, J);
    VISCOUS viscous(N, Re, xvel, yvel, zx, ey, J);
    for (int k = 0; k < 3; ++k)
    {
        for (int i = 1; i < N - 1; ++i)
        {
            for (int j = 1; j < N - 1; ++j)
            {
                rhs[k][i][j] = 0.5 * J[i][j] * ((convective.Estar1[k][i][j + 1] - convective.Estar1[k][i][j - 1]) + (convective.Estar2[k][i + 1][j] - convective.Estar2[k][i - 1][j]))
                               - (1 / Re) * J[i][j] * ((viscous.Estarv1[k][i][j + 1] - viscous.Estarv1[k][i][j - 1]) + (viscous.Estarv2[k][i + 1][j] - viscous.Estarv2[k][i - 1][j]));
            }
        }
    }
    spectralradius(N, J, convective.U, convective.V, viscous.g11, viscous.g22);
    Diss( N, Press, xvel, yvel);
    smoothing(N);
}


void RHS::spectralradius(   const int& N,
                            const vector<vector<double> > &J,
                            const vector<vector<double> > &U,
                            const vector<vector<double> > &V,
                            const vector<vector<double> > &g11,
                            const vector<vector<double> > &g22)
{
    Resize2D(N, rho1); Resize2D(N, rho2);

    for (int i = 0; i < N; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++)
        {
            rho1[i][j] = (1 / J[i][j]) * (abs(U[i][j]) + sqrt((pow(U[i][j], 2) + g11[i][j])));
            rho2[i][j] = (1 / J[i][j]) * (abs(V[i][j]) + sqrt((pow(V[i][j], 2) + g22[i][j])));
        }
    }
}

void RHS::Diss( const int& N,
                const vector<vector<double> > &Press,
                const vector<vector<double> > &xvel,
                const vector<vector<double> > &yvel)
{
    Resize2D(N, D1_c); Resize2D(N, D2_c); Resize2D(N, D1_u); Resize2D(N, D2_u); Resize2D(N, D1_v); Resize2D(N, D2_v);
    Resize2D(N, Diss_c); Resize2D(N, Diss_u); Resize2D(N, Diss_v);
    eps = 0.01;
    for (int i = 2; i < N - 2; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 2; j < N - 2; j++)
        {
            D1_c[i][j] = eps * (rho1[i][j] + rho1[i + 1][j]) * 0.5 * (Press[i + 2][j] - 3 * Press[i + 1][j] + 3 * Press[i][j] - Press[i - 1][j]);
            D1_u[i][j] = eps * (rho1[i][j] + rho1[i + 1][j]) * 0.5 * (xvel[i + 2][j] - 3 * xvel[i + 1][j] + 3 * xvel[i][j] - xvel[i - 1][j]);
            D1_v[i][j] = eps * (rho1[i][j] + rho1[i + 1][j]) * 0.5 * (yvel[i + 2][j] - 3 * yvel[i + 1][j] + 3 * yvel[i][j] - yvel[i - 1][j]);
            D2_c[i][j] = eps * (rho2[i][j] + rho2[i][j + 1]) * 0.5 * (Press[i][j + 2] - 3 * Press[i][j + 1] + 3 * Press[i][j] - Press[i][j - 1]);
            D2_u[i][j] = eps * (rho2[i][j] + rho2[i][j + 1]) * 0.5 * (xvel[i][j + 2] - 3 * xvel[i][j + 1] + 3 * xvel[i][j] - xvel[i][j - 1]);
            D2_v[i][j] = eps * (rho2[i][j] + rho2[i][j + 1]) * 0.5 * (yvel[i][j + 2] - 3 * yvel[i][j + 1] + 3 * yvel[i][j] - yvel[i][j - 1]);
            Diss_c[i][j] = ((D1_c[i + 1][j] + D1_c[i][j]) * 0.5 - 0.5 * (D1_c[i][j] + D1_c[i - 1][j])) + ((D2_c[i + 1][j] + D2_c[i][j]) * 0.5 - 0.5 * (D2_c[i][j] + D2_c[i][j - 1]));
            Diss_u[i][j] = ((D1_u[i + 1][j] + D1_u[i][j]) * 0.5 - 0.5 * (D1_u[i][j] + D1_u[i - 1][j])) + ((D2_u[i + 1][j] + D2_u[i][j]) * 0.5 - 0.5 * (D2_u[i][j] + D2_u[i][j - 1]));
            Diss_v[i][j] = ((D1_v[i + 1][j] + D1_v[i][j]) * 0.5 - 0.5 * (D1_v[i][j] + D1_v[i - 1][j])) + ((D2_v[i + 1][j] + D2_v[i][j]) * 0.5 - 0.5 * (D2_v[i][j] + D2_v[i][j - 1]));

        }
    }
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; i++)
    {
        Diss_c[i][0] = Diss_c[i][1];
        Diss_u[i][0] = Diss_u[i][1];
        Diss_v[i][0] = Diss_v[i][1];
        Diss_c[i][N - 2] = Diss_c[i][N - 3];
        Diss_u[i][N - 2] = Diss_u[i][N - 3];
        Diss_v[i][N - 2] = Diss_v[i][N - 3];
        Diss_c[i][N - 1] = Diss_c[i][N - 2];
        Diss_u[i][N - 1] = Diss_u[i][N - 2];
        Diss_v[i][N - 1] = Diss_v[i][N - 2];
    }
    #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < N; j++)
    {
        Diss_c[0][j] = (Diss_c[1][j]);
        Diss_u[0][j] = (Diss_u[1][j]);
        Diss_v[0][j] = (Diss_v[1][j]);
        Diss_c[N - 2][j] = Diss_c[N - 3][j];
        Diss_u[N - 2][j] = Diss_u[N - 3][j];
        Diss_v[N - 2][j] = Diss_v[N - 3][j];
        Diss_c[N - 1][j] = Diss_c[N - 2][j];
        Diss_u[N - 1][j] = Diss_u[N - 2][j];
        Diss_v[N - 1][j] = Diss_v[N - 2][j];
    }
    for (int i = 0; i < N; i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; j++)
        {
            rhs[0][i][j] = rhs[0][i][j] + Diss_c[i][j];
            rhs[1][i][j] = rhs[1][i][j] + Diss_u[i][j];
            rhs[2][i][j] = rhs[2][i][j] + Diss_v[i][j];
        }
    }

}

void RHS::smoothing(const int& N)
{
    ep = 0.01;
    vector<vector<double> > Dummyu (N, vector<double>(N, 0));
    vector<vector<double> > Dummyv (N, vector<double>(N, 0));
    vector<vector<double> > Dummyc (N, vector<double>(N, 0));

    vector<double> a(N);
    vector<double> b(N);
    vector<double> c(N);
    vector<double> du(N);
    vector<double> dv(N);
    vector<double> dc(N);


    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; ++i)
    {
        a[i] = 0;
        b[i] = 0;
        c[i] = 0;
        dc[i] = 0;
        du[i] = 0;
        dv[i] = 0;
        for (int j = 0; j < N; ++j)
        {
            Dummyc[i][j] = 0;
            Dummyu[i][j] = 0;
            Dummyv[i][j] = 0;
        }
    }

    a[0] = 0;
    b[0] = 1 + 2 * ep;
    c[0] = -ep;
    a[N - 1] = -ep;
    b[N - 1] = 1 + 2 * ep;
    c[N - 1] = 0;
    for (int i = 0; i < N; ++i)
    {
        dc[0] = rhs[0][i][0];
        du[0] = rhs[1][i][0];
        dv[0] = rhs[2][i][0];
        #pragma omp parallel for schedule(dynamic)
        for (int j = 1; j < N - 1; ++j)
        {
            a[j] = -ep;
            b[j] = 1 + 2 * ep;
            c[j] = -ep;
            dc[j] = rhs[0][i][j];
            du[j] = rhs[1][i][j];
            dv[j] = rhs[2][i][j];
        }
        dc[N - 1] = rhs[0][i][N - 1];
        du[N - 1] = rhs[1][i][N - 1];
        dv[N - 1] = rhs[2][i][N - 1];
        TRI(1, N - 1, a, b, c, dc);
        TRI(1, N - 1, a, b, c, du);
        TRI(1, N - 1, a, b, c, dv);
        dv[0] = rhs[2][i][0];
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; ++j)
        {
            Dummyc[i][j] = dc[j];
            Dummyu[i][j] = du[j];
            Dummyv[i][j] = dv[j];
        }
    }

    for (int i = 0; i < N; ++i)
    {
        dc[0] = Dummyc[i][0];
        du[0] = Dummyu[i][0];
        dv[0] = Dummyv[i][0];
        #pragma omp parallel for schedule(dynamic)
        for (int j = 1; j < N - 1; ++j)
        {
            a[j] = -ep;
            b[j] = 1 + 2 * ep;
            c[j] = -ep;
            dc[j] = Dummyc[i][j];
            du[j] = Dummyu[i][j];
            dv[j] = Dummyv[i][j];
        }
        dc[N - 1] = Dummyc[i][N - 1];
        du[N - 1] = Dummyc[i][N - 1];
        dv[N - 1] = Dummyc[i][N - 1];
        TRI(1, N - 1, a, b, c, dc);
        TRI(1, N - 1, a, b, c, du);
        TRI(1, N - 1, a, b, c, dv);
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < N; ++j)
        {
            rhs[0][i][j] = dc[j];
            rhs[1][i][j] = du[j];
            rhs[2][i][j] = dv[j];
        }
    }
}
