/*
this is the c++ code for the DNO in Navier's equation. 
*/
#include<iostream>
#include<complex>
#include<math.h>

#define pi 4.0*atan(1)
using namespace std;

void dg_3d_scatter_pre(double alpha,double beta,double* alphap,double* betap,double* pp1,double* pp2)
{

	double pp1pre[nx1],pp2pre[nx2];
	

}

int main()
{
/*
% N - Number of Taylor orders
% nx1 - Number of equally spaced gridpoints on [0,d1]
% nx2 - Number of equally spaced gridpoints on [0,d2]
*/
	int N,nx1,nx2;
	double epsilon,d1,d2,dx1,dx2;
	N = 2;
	nx1 = 8;
	nx2 = 8;
	epsilon =1e-6;
	d1 = 0.53;
	d2 = d1;

	dx1 = d1/nx1;
	dx2 = d2/nx2;

	double *xx1,*xx2;
	xx1 = (double*)calloc(nx1*nx2,sizeof(double));
	xx2 = (double*)calloc(nx1*nx2,sizeof(double));
	for(int i=0;i<nx1;i++)
		for(int l=0;l<nx2;l++)
		{
			xx1[nx2*i + l] = 0;
			xx2[nx2*i + l] = 0;
		}

	for(int i=0;i<nx1;i++)
		for(int l=0;l<nx2;l++)
		{
			xx1[nx2*i + l] = (i-1)*dx1;
			xx2[nx2*i + l] = (l-1)*dx2;
		}
/*
	% mu lambda rho- physical parameters
	% c_1, c_2 are square of the speed of dilatational and rotational wave
	% omega - is harmonic period of time
*/
	double mu,lambda,rho,omega,alpha1,alpha2,c_1,c_2,k_1,k_2,beta1,beta2;
	mu = 1;
	lambda = 1;
	rho = 1;
	omega = 14;
	alpha1 = 0.0;
	alpha2 = alpha1;
	c_1 = (lambda + 2*mu)/rho;
	c_2 = mu/rho;
	k_1 = sqrt(omega/c_1);
	k_2 = sqrt(omega/c_2);
	beta1 = k_1;
	beta2 = k_2;				

// 3D scatter
	double *alpha1p,*alpha2p,*pp1,*pp2;
	complex<double> *betap1,*betap2;
	alpha1p = (double*)calloc(nx1*nx2,sizeof(double));
	alpha2p = (double*)calloc(nx1*nx2,sizeof(double));
	pp1 = (double*)calloc(nx1*nx2,sizeof(double));
	pp2 = (double*)calloc(nx1*nx2,sizeof(double));
	for(int i=0;i<nx1;i++)
		for(int l=0;l<nx2;l++)
		{
			pp1[nx2*i + l] = 0;
			pp2[nx2*i + l] = 0;
			alpha1p[nx2*i + l] = 0;
			alpha2p[nx2*i + l] = 0;
			betap1[nx2*i + l] = 0;
			betap2[nx2*i + l] = 0;
		}


	dg_3d_scatter_pre(alpha1,alpha2,alpha1p,alpha2p,pp1,pp2);






	 

}

