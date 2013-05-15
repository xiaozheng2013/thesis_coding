/*
this is the c++ code for the DNO in Navier's equation. 
*/
#include<iostream>
#include<fstream>
#include<iomanip>
#include<complex>
#include<fftw3.h>
#include<math.h>

#define pi 4.0*atan(1)
using namespace std;
complex<double> one(0,1);

void init(complex<double>* a, int b)
{
	for(int i=0;i<b;i++)
		a[i] = 0.0;
}

void fc2c(fftw_complex *in, complex<double> *out, int n)  //http://www.physicsforums.com/showthread.php?t=97608
{
	int i;
	for(i=0;i<n;i++) 
	{
		out[i] = complex<double>(in[i][0],in[i][1]);
	}
}

void c2fc(complex<double> *in, fftw_complex *out, int n)
{
	int i;
	for(i=0;i<n;i++) 
	{
		out[i][0] = real(in[i]);
		out[i][1] = imag(in[i]); 
	}
}

void fft2(complex<double> *f,int nx,int ny,complex<double> *fhat)
{
	fftw_complex *in,*out;
	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
	 c2fc(f,in,nx*ny);

	fftw_plan p;
	p = fftw_plan_dft_2d(nx,ny,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	
	fc2c(out,fhat,nx*ny);
	fftw_destroy_plan(p);
}

void ifft2(complex<double> *fhat,int nx,int ny,complex<double> *f)
{
	fftw_complex *in,*out;
	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(nx*ny));
	 c2fc(fhat,in,nx*ny);

	fftw_plan p;
	p = fftw_plan_dft_2d(nx,ny,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	
	fc2c(out,f,nx*ny);
	fftw_destroy_plan(p);
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
		{
			int k= i*ny + j;
			f[k] = f[k]/(double)(nx*ny);
		}
}

void dg_3d_scatter_pre(int nx1,int nx2, double d1,double d2, double alpha,double beta,double* alphap,double* betap,double* pp1,double* pp2)
{

	double pp1pre[nx1],pp2pre[nx2];
	for(int i=0;i<nx1;i++)
	{
		if(i<nx1/2)
		pp1pre[i] = (2*pi)/d1*i;
		else
		pp1pre[i] = (2*pi)/d1*(i-nx1); 
	}
	for(int i=0;i<nx2;i++)
	{
		if(i<nx2/2)
		pp2pre[i] = (2*pi)/d2*i;
		else
		pp2pre[i] = (2*pi)/d2*(i-nx2); 
	}

	for(int i=0;i<nx1;i++)
		for(int j=0;j<nx2;j++)
		{
			int k = i*nx2 + j; 
			pp1[k] = pp1pre[i];
			pp2[k] = pp2pre[j];
			alphap[k] = alpha + pp1[k];
			betap[k] = beta + pp2[k];
		}

}

void dg_3d_scatter(int nx1,int nx2,double alpha1,double alpha2,double beta1,double* alphap,double* betap,complex<double>* betap1)
{
	double K = sqrt(pow(alpha1,2) + pow(alpha2,2) + pow(beta1,2));	
	for(int i=0;i<nx1;i++)
		for(int j=0;j<nx2;j++)
		{
			int k = i*nx2 + j;
			if(pow(alphap[k],2) + pow(betap[k],2) < pow(K,2))
				betap1[k] = sqrt(pow(K,2) - pow(alphap[k],2) - pow(betap[k],2));
			else
				betap1[k] = one*sqrt(pow(alphap[k],2) + pow(betap[k],2) - pow(K,2));
		}
}


void fhatsetup3d(int nx,int ny,complex<double>* fhat)
{
	switch(1)
		case 1:
		{
			int k1 = ny + 1;
			fhat[k1] = nx*ny/2.0;
			int k2 = (nx-1)*ny + (ny-1);
			fhat[k2] = nx*ny/2.0;
			break;

		}
}

void newconv3d(int nx,int ny,complex<double>* a, complex<double>* b, complex<double>* c)
{
	complex<double> *a_tmp,*b_tmp;	
	a_tmp = (complex<double>*)calloc(nx*ny,sizeof(complex<double>));
	b_tmp = (complex<double>*)calloc(nx*ny,sizeof(complex<double>));
	ifft2(a,nx,ny,a_tmp);
	ifft2(b,nx,ny,b_tmp);
	for(int i=0;i<nx*ny;i++)
	{
		a_tmp[i] *= b_tmp[i];
	}
	fft2(a_tmp,nx,ny,c);
}

void cpnsetup3d(int nx,int ny,complex<double>* fhat,int N,complex<double>* fpn)
{
	fpn[0] = nx*ny;
	for(int i=0;i<nx*ny;i++)
	{
		fpn[i+nx*ny] = fhat[i];	
	}
	for(int i=1;i<N;i++)
	{
		int ind = i*nx*ny;
		int ind1 = (i+1)*nx*ny;
		newconv3d(nx,ny,fhat,&fpn[ind],&fpn[ind1]);
		for(int j=0;j<nx*ny;j++)
			fpn[ind1+j] /= (double)(i+1);
	}

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
	betap1 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	betap2 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));

	for(int i=0;i<nx1;i++)
		for(int l=0;l<nx2;l++)
		{
			int k = nx2*i + l;
			pp1[k] = 0;
			pp2[k] = 0;
			alpha1p[k] = 0;
			alpha2p[k] = 0;
			betap1[k] = 0;
			betap2[k] = 0;
		}

	dg_3d_scatter_pre(nx1,nx2,d1,d2,alpha1,alpha2,alpha1p,alpha2p,pp1,pp2);
	dg_3d_scatter(nx1,nx2,alpha1,alpha2,beta1,alpha1p,alpha2p,betap1);
	dg_3d_scatter(nx1,nx2,alpha1,alpha2,beta2,alpha1p,alpha2p,betap2);
// finished checking at the 1st checking point 5/12/13

// set up the powers of f
	complex<double> *fpn,*fhat,*f,*fx1hat,*fx1,*fx2hat,*fx2;
	fpn = (complex<double> *)calloc(nx1*nx2*max(2,N+1),sizeof(complex<double>));
	fhat = (complex<double> *)calloc(nx1*nx2,sizeof(complex<double>));
	f = (complex<double> *)calloc(nx1*nx2,sizeof(complex<double>));
	fx1hat = (complex<double> *)calloc(nx1*nx2,sizeof(complex<double>));
	fx1 = (complex<double> *)calloc(nx1*nx2,sizeof(complex<double>));
	fx2hat = (complex<double> *)calloc(nx1*nx2,sizeof(complex<double>));
	fx2 = (complex<double> *)calloc(nx1*nx2,sizeof(complex<double>));
	for(int i=0;i<nx1;i++)
		for(int j=0;j<nx2;j++)
		{
			int k = i*nx2 + j;
			fhat[k] = 0.0;			 
			f[k] = 0.0;		
			fx1hat[k] = 0.0;		
			fx1[k] = 0.0;		
			fx2hat[k] = 0.0;		
			fx2[k] = 0.0;		
		}
	for(int i=0;i<max(2,N+1);i++)
		for(int j=0;j<nx1;j++)
			for(int l=0;l<nx2;l++)
			{
				int k = i*nx1*nx2 + j*nx2 + l;
				fpn[k] = 0.0;
			}
	fhatsetup3d(nx1,nx2,fhat);
	ifft2(fhat,nx1,nx2,f);
	for(int i=0;i<nx1;i++)
		for(int j=0;j<nx2;j++)
		{
			int k = i*nx2 + j;
			fx1hat[k] = one*pp1[k]*fhat[k];
			fx2hat[k] = one*pp2[k]*fhat[k];
		}
	ifft2(fx1hat,nx1,nx2,fx1);
	ifft2(fx2hat,nx1,nx2,fx2);
	
	cpnsetup3d(nx1,nx2,fhat,N,fpn);
	
// checkpoint 2, finished checking on 2013/5/14

// form zeta(Dirichlet) and psi(Neumann) data
	complex<double> *N_inv,*N_orig,*a1,*a2,*a3;
	N_inv = (complex<double>*)calloc(3*3*nx1*nx2,sizeof(complex<double>));
	N_orig = (complex<double>*)calloc(3*3*nx1*nx2,sizeof(complex<double>));
	a1 = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
	a2 = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
	a3 = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
	init(N_inv,3*3*nx1*nx2);
	init(N_orig,3*3*nx1*nx2);
	init(a1,3*nx1*nx2);
	init(a2,3*nx1*nx2);
	init(a3,3*nx1*nx2);
	

	ofstream o1("fpn.dat");
	for(int l=0;l<max(2,N+1);l++)	
	{
	for(int i=0;i<nx1;i++)
	{
		for(int j=0;j<nx2;j++)
		{ 
			int k=l*nx1*nx2 + i*nx2 + j;
			o1<<setprecision(20)<<fpn[k]<<" ";
		}
		o1<<endl;
	}
	o1<<endl;
	}




	o1.close();	 
	delete pp1;
	delete pp2;
	delete alpha1p;
	delete alpha2p;
	delete betap1;
	delete betap2;

	delete fhat;
	delete fx1;
	delete fx2;
	delete fx2hat;
	delete fx1hat;
	delete fpn;
	return 0;

}

