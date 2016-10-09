

void BC (int N, double** u_new,double** v_new,double** p_new,double** dummyu)
{	
	for (int i=0; i<N; i++)
	{
		u_new[i][0]=dummyu[i][0];
		u_new[0][i]=0;
		u_new[N-1][i]=0;
		u_new[i][N-1]=u_new[i][N-2];

		v_new[i][0]=0;
		v_new[0][i]=0;
		v_new[N-1][i]=0;
		v_new[i][N-1]=v_new[i][N-2];

		p_new[0][i]=p_new[1][i];
		p_new[i][0]=p_new[i][1];
		p_new[N-1][i]=p_new[N-2][i];
		p_new[i][N-1]=p_new[i][N-2];
	}
	
}
