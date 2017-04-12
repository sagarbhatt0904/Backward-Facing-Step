/*

	Class to initialize channel flow. Currently it supperts 2D normalized channel flow. Other functionality to be added in future.
	Created on: April 5, 2017
    Author: Sagar Bhatt

*/


#include <vector>
#include <string>

using namespace std;

class INIT
{

public:
	INIT() {}
	INIT(int& N, vector<double> &x, vector<double> &y, vector<vector<double> > &xvel, vector<vector<double> > &yvel, vector<vector<double> > &Press);
};



INIT::INIT(int& N, vector<double> &x, vector<double> &y, vector<vector<double> > &xvel, vector<vector<double> > &yvel, vector<vector<double> > &Press)
{

	/* Initializing U,V, P*/
	int kl = 0;
	vector<vector<double> > xd (N, vector<double>(N, 0));

	for (int j = ((N - 1) / 2); j < N ; j++)
	{
		xd[1][j] = 1 - ((y[j] - 1) * (y[j] - 1));
	}
	for (int j = 0; j < (N - 1) / 2; j++)
	{
		xd[1][j] = xd[1][(N - 1) - kl];
		kl = kl + 1;
	}

	xd[1][(N - 1) / 2] = 1;
	for (int j = (N - 1) / 2; j < N; j++)
	{
		xvel[j][1] = xd[1][(2 * j) - (N - 1)];
	}
}