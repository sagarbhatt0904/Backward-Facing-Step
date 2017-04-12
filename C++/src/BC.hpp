
using namespace std;

void BC (const int& N, vector<vector<double> > &u_new,vector<vector<double> > &v_new,vector<vector<double> > &p_new)
{	
	for (int i=0; i<N; i++)
	{
		u_new[i][0]=0.5*(u_new[i][1]+u_new[i][N-2]);
		u_new[0][i]=0.5*(u_new[1][i]+u_new[N-2][i]);
		u_new[N-1][i]=u_new[0][i];
		u_new[i][N-1]=u_new[i][0];

		v_new[i][0]=0.5*(v_new[i][1]+v_new[i][N-2]);
		v_new[0][i]=0.5*(v_new[1][i]+v_new[N-2][i]);
		v_new[N-1][i]=v_new[0][i];
		v_new[i][N-1]=v_new[i][0];

		p_new[i][0]=0.5*(p_new[i][1]+p_new[i][N-2]);
		p_new[0][i]=0.5*(p_new[1][i]+p_new[N-2][i]);
		p_new[N-1][i]=p_new[0][i];
		p_new[i][N-1]=p_new[i][0];
	}
	
}
