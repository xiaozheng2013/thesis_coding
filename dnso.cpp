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

void mult(int nx,complex<double> a, complex<double> *x, complex<double> *y)
{
	for(int i=0;i<nx;i++)
		y[i] = a*x[i];
}

void mult(int nx,double* a, complex<double> *x, complex<double> *y)
{
	for(int i=0;i<nx;i++)
		y[i] = a[i]*x[i];
}

double radius ( complex<double> c )
// Returns the radius of the complex number.
{
   double result,re,im;
   re = real(c);im=imag(c);
   result = re*re + im*im;
   return sqrt(result);
}

void init(complex<double>* a, int b)
{
	for(int i=0;i<b;i++)
		a[i] = 0.0;
}

void init(double* a, int b)
{
	for(int i=0;i<b;i++)
		a[i] = 0.0;
}

void fc2c(fftw_complex *in, complex<double> *out, int n)  //http://www.physicsforums.com/showthread.php?t=97608
{
	int i;
	for(i=0;i<n;i++) {
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

void complexBackSubstitution(complex<double> *U[],
                     complex<double> *b,
                     complex<double> *x,
                     int M)
{
  int j,k,n;
  complex<double> sum;


  for(j=M-1;j>=0;j--){
    sum = 0.0;
    for(k=j+1;k<M;k++){
      sum += U[j][k]*x[k];
      }
    x[j]=(b[j]-sum)/U[j][j];
    }

} 

void complexTriangulate(complex<double> *A[],
                complex<double> *b,
                complex<double> *C[],
                complex<double> *d,
                int M)
{
  complex<double> temprow[M];
  complex<double> potentialpivot,pivot,factor,tempd;
  int j,k,m,n,pivotrow;

  //n = b.Length();
  n = M;
  for(j=0;j<M;j++)
  {
      temprow[j] = 0.0;
      d[j] = b[j];
      for(m=0;m<M;m++)
          C[j][m] = A[j][m];
  }
 	for(j=0;j<M-1;j++){

    /* Find largest pivot */
    pivot = C[j][j];
    pivotrow = j;
    for(m=j+1;m<n;m++){
      potentialpivot = C[m][j];
      //      if(fabs(potentialpivot)>fabs(pivot)){
      double judge;
      judge = radius(potentialpivot) - radius(pivot);
      if( judge > 1e-10 ){
    pivot = potentialpivot;
    pivotrow = m;
    }
      }

    /* Interchange rows and RHS */
    for(m=0;m<n;m++){
      temprow[m] = C[j][m];
      C[j][m] = C[pivotrow][m];
      C[pivotrow][m] = temprow[m];
      }
 	tempd = d[j];
    d[j] = d[pivotrow];
    d[pivotrow] = tempd;

    /* Resume factorization */
    for(k=j+1;k<M;k++){
      factor = C[k][j]/C[j][j];
      for(m=j;m<n;m++)
      {

          C[k][m] = C[k][m] - factor * C[j][m];
      }
      d[k] = d[k] - factor * d[j];
      }
    }

} 

void complexGaussElimination(complex<double> *A[],
                     complex<double> *b,
                     complex<double> *x, int M)
{
  complex<double> *C[M];
  complex<double> d[M];

  for(int i=0;i<M;i++)
  {
    d[i] = 0.0;
    C[i] = (complex<double> *)calloc(M,sizeof(complex<double>));
    for(int j = 0;j<M;j++)
        C[i][j] = 0.0;
  }

  complexTriangulate(A,b,C,d,M);
  complexBackSubstitution(C,d,x,M);

}  
			
void stress_0( complex<double>* Un, double alpha1p, double alpha2p, complex<double> betap1,complex<double> betap2, double lambda, double mu, int nx1, int nx2, complex<double>* right)
{
	complex<double> *wave,*p_wave,*s_wave;
	wave = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	p_wave = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	s_wave = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	init(wave,3*nx1*nx2);
	init(p_wave,3*nx1*nx2);
	init(s_wave,3*nx1*nx2);

	mult(3*nx1*nx2,1,Un,p_wave);
	mult(3*nx1*nx2,1,&Un[3*nx1*nx2],s_wave);
	mult(3*nx1*nx2,1,&Un[6*nx1*nx2],wave);
	for(int i=0;i<nx1*nx2*3;i++)
	{
		s_wave[i] += wave[i];
		wave[i] = p_wave[i] + s_wave[i];
	}

	// useful patterns, which are all in F-S
	complex<double> *DIV,*par_1_3,*par_2_3,*par_3_3;
	DIV = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	par_1_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	par_2_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	par_3_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	init(DIV,nx1*nx2);
	init(par_1_3,nx1*nx2);
	init(par_2_3,nx1*nx2);
	init(par_3_3,nx1*nx2);
	
	for(int i=0;i<nx1;i++)
		for(int j=0;j<nx2;j++)
		{
			int k = i*nx2 + j;
			int ind2 = nx1*nx2;
			int ind3 = 2*nx1*nx2;

			DIV[k] = one*(alpha1p*p_wave[k] + alpha2p*p_wave[ind2 + k] + betap1*p_wave[ind3 + k]);

			par_1_3[k] = one*(alpha1p*wave[ind3 + k] + betap1*p_wave[k] + betap2*s_wave[k]);
			par_2_3[k] = one*(alpha2p*wave[ind3 + k] + betap1*p_wave[ind2 + k] + betap2*s_wave[ind2 + k]);
			par_3_3[k] = one*(betap1*p_wave[ind3 + k] + betap2*s_wave[ind3 + k]);

		}	
		
	int ind2 = nx1*nx2;
	mult(nx1*nx2,-mu,par_1_3,right);
	mult(nx1*nx2,-mu,par_2_3,&right[ind2]);
	for(int i=0;i<nx1;i++)
		for(int j=0;j<nx2;j++)
		{
			int k = i*nx2 + j;
			int ind = 2*nx1*nx2;
			right[ind + k] = -lambda*DIV[k] - 2*mu*par_3_3[k];

		}	

} 
	
void g0_engine(int n,int nx1,int nx2,complex<double>* Un1,complex<double>* Un2,complex<double>* Un3,complex<double>* fpn,double* alpha1p,double* alpha2p,complex<double>* betap1,complex<double>* betap2,double lambda,double mu,complex<double>* G0_I_J1,complex<double>* G0_I_J2,complex<double>* G0_I_J3) 
{

	for(int i=0;i<nx1;i++)
		for(int j=0;j<nx2;j++)
		{
			complex<double> *right1,*right2,*right3;
			right1 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
			right3 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
			right2 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
			init(right1,3*nx1*nx2);
			init(right2,3*nx1*nx2);
			init(right3,3*nx1*nx2);

			int ind = nx1*nx2*(i*nx2*9 + j*9);
			int ind1 = i*nx2 + j;

			stress_0( &Un1[ind], alpha1p[ind1],alpha2p[ind1],betap1[ind1],betap2[ind1],lambda,mu,nx1,nx2, right1); 
			stress_0( &Un2[ind], alpha1p[ind1],alpha2p[ind1],betap1[ind1],betap2[ind1],lambda,mu,nx1,nx2, right2); 
			stress_0( &Un3[ind], alpha1p[ind1],alpha2p[ind1],betap1[ind1],betap2[ind1],lambda,mu,nx1,nx2, right3); 

			int ind2 = nx1*nx2*(i*nx2*3 + j*3);
			mult(3*nx1*nx2,1,right1,&G0_I_J1[ind2]);
			mult(3*nx1*nx2,1,right2,&G0_I_J2[ind2]);
			mult(3*nx1*nx2,1,right3,&G0_I_J3[ind2]);

			delete right1,right2,right3;
		}	

}

void stress_pre(int nx1,int nx2,complex<double>* Un,double alpha1p,double alpha2p,complex<double> betap1,complex<double> betap2,complex<double>* p_wave,complex<double>* s_wave,complex<double>* wave,complex<double>* div,complex<double>* par12,complex<double>* par13,complex<double>* par23,complex<double>* par33)
{
/*
% this code did all the preparation work for calculation of DNO at each
% order, get p-wave,s-wave, and many useful patterns before reall
% calculation.
*/
	int ind2 = 3*nx1*nx2,ind3 = 6*nx1*nx2;
	mult(3*nx1*nx2,1,Un,p_wave);
	mult(3*nx1*nx2,1,&Un[ind2],s_wave);
	mult(3*nx1*nx2,1,&Un[ind3],wave);
	for(int i=0;i<3*nx2*nx1;i++)
	{
		s_wave[i] += wave[i];
		wave[i] = p_wave[i] + s_wave[i];
	}

	ind2 = nx1*nx2;
	ind3 = 2*nx2*nx1;
	for(int i=0;i<nx1*nx2;i++)	
	{
		
		div[i] = one*alpha1p*p_wave[i] + one*alpha2p*p_wave[ind2 + i] 
				+ one*betap1*p_wave[ind2 + i];
		par12[i] = one*(alpha1p*wave[ind2 + i] + alpha2p*wave[i]);
		par13[i] = one*(alpha1p*wave[ind3 + i] + betap1*p_wave[i] 
					+ betap2*s_wave[i]);
		par23[i] = one*(alpha2p*wave[ind3 + i] + betap1*p_wave[ind2 + i]
					+ betap2*s_wave[ind2 + i]);
		par33[i] = one*(betap1*p_wave[ind3 + i] + betap2*s_wave[ind3 + i]);
	}

}

void stress_n(complex<double> *u1,complex<double> *u2,double alpha1p,double alpha2p,complex<double> betap1,complex<double> betap2,double lambda,double mu,int nx1,int nx2,complex<double> *fx,complex<double> *fy,complex<double> *right)
{
/*
% This function calculates the n-th order perturbation FC of the stress,
% which is on the right hand side of DNO recursive formula
*/
	complex<double> *p_wave_1,*s_wave_1,*wave_1;
	complex<double> *p_wave_2,*s_wave_2,*wave_2;
	p_wave_1 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	s_wave_1 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	wave_1 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	p_wave_2 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	s_wave_2 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	wave_2 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));

	init(p_wave_1,3*nx1*nx2);		
	init(s_wave_1,3*nx1*nx2);		
	init(wave_1,3*nx1*nx2);		
	init(p_wave_2,3*nx1*nx2);		
	init(s_wave_2,3*nx1*nx2);		
	init(wave_2,3*nx1*nx2);		

	complex<double> *DIV_1,*par_1_2,*par_1_3,*par_2_3,*par_3_3;
	DIV_1 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	par_1_2 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	par_1_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	par_2_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	par_3_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	init(DIV_1,nx1*nx2);		
	init(par_1_2,nx1*nx2);		
	init(par_1_3,nx1*nx2);		
	init(par_2_3,nx1*nx2);		
	init(par_3_3,nx1*nx2);		

	complex<double> *DIV_2,*Par_1_2,*Par_1_3,*Par_2_3,*Par_3_3;
	DIV_2 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	Par_1_2 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	Par_1_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	Par_2_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	Par_3_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	init(DIV_2,nx1*nx2);		
	init(Par_1_2,nx1*nx2);		
	init(Par_1_3,nx1*nx2);		
	init(Par_2_3,nx1*nx2);		
	init(Par_3_3,nx1*nx2);		

	stress_pre(nx1,nx2,u1,alpha1p,alpha2p,betap1,betap2,p_wave_1,s_wave_1,wave_1,DIV_1,par_1_2,par_1_3,par_2_3,par_3_3);
	stress_pre(nx1,nx2,u2,alpha1p,alpha2p,betap1,betap2,p_wave_2,s_wave_2,wave_2,DIV_2,Par_1_2,Par_1_3,Par_2_3,Par_3_3);
		
	complex<double> *tmp1,*tmp2,*tmp3;
	tmp1 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	tmp2 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	tmp3 = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
	init(tmp1,nx1*nx2*3);
	init(tmp2,nx1*nx2*3);
	init(tmp3,nx1*nx2*3);

// now calculate the first dimension
	mult(nx1*nx2,lambda,DIV_1,tmp1);
	mult(nx1*nx2,2*mu*one*alpha1p,wave_1,tmp2);
	for(int i=0;i<nx1*nx2;i++)
	{
		tmp1[i] += tmp2[i];
	}
	
	mult(nx1*nx2,mu,par_1_2,tmp2);
	mult(nx1*nx2,mu,Par_1_2,tmp3);

	newconv3d(nx1,nx2,tmp1,fx,tmp1);
	newconv3d(nx1,nx2,tmp2,fy,tmp2);
	for(int i=0;i<nx1*nx2;i++)
	{
		right[i] = tmp1[i] + tmp2[i] - tmp3[i];
	} 

// now calculate the second dimension
	int ind2 = nx1*nx2;
	mult(nx1*nx2,mu,par_1_2,&tmp1[ind2]);
	mult(nx1*nx2,lambda,DIV_1,&tmp2[ind2]);
	mult(nx1*nx2,2*mu*one*alpha2p,&wave_1[ind2],&tmp3[ind2]);
	for(int i=0;i<nx1*nx2;i++)
	{
		tmp2[ind2 + i] += tmp3[ind2 + i];
	}
	mult(nx1*nx2,mu,Par_2_3,&tmp3[ind2]);

	newconv3d(nx1,nx2,&tmp1[ind2],fx,&tmp1[ind2]);
	newconv3d(nx1,nx2,&tmp2[ind2],fy,&tmp2[ind2]);
	for(int i=ind2;i<ind2 + nx1*nx2;i++)
	{
		right[i] = tmp1[i] + tmp2[i] - tmp3[i];
	} 
// now calculate the 3rd dimension
	int ind3 = 2*nx2*nx1;
	mult(nx1*nx2,mu,par_1_3,&tmp1[ind3]);
	mult(nx1*nx2,lambda,DIV_2,&tmp2[ind3]);
	mult(nx1*nx2,2*mu,Par_3_3,&tmp3[ind3]);
	for(int i=ind3;i<ind3 + nx1*nx2;i++)
	{
		tmp3[i] += tmp2[i];
	}
	mult(nx1*nx2,mu,par_2_3,&tmp2[ind3]);
	
	newconv3d(nx1,nx2,&tmp1[ind3],fx,&tmp1[ind3]);
	newconv3d(nx1,nx2,&tmp2[ind3],fy,&tmp2[ind3]);
	for(int i=ind3;i<ind3 + nx1*nx2;i++)
	{
		right[i] = tmp1[i] + tmp2[i] - tmp3[i];
	} 

	delete p_wave_1,s_wave_1,wave_1;
	delete p_wave_2,s_wave_2,wave_2;
	delete DIV_1,par_1_2,par_1_3,par_2_3,par_3_3;
	delete DIV_2,Par_1_2,Par_1_3,Par_2_3,Par_3_3;
	delete tmp1,tmp2,tmp3;

}

void gn_engine(int n,int nx1,int nx2,complex<double>* Un1,complex<double>* Un2,complex<double>* Un3,complex<double>* fpn,double* alpha1p,double* alpha2p,complex<double>* betap1,complex<double>* betap2,double lambda,double mu,complex<double>* Gn_m_mUm1,complex<double>* Gn_m_mUm2,complex<double>* Gn_m_mUm3) 
{
/*	% 
% this function will solve all the G_n[U(i,j)]
% 
% Inputs:
% . n - the order of DNSO
% . Un - Fourier coefficients of Dirichlet data of previous n orders
% . fpn - Fourier coefficients of f^n/n!
% . alphap - quasiwavenumbers: \alpha + (2 \pi/d) p
% . betap - \beta_p = \sqrt{ k^2 - \alpha_p^2 }
% . Gn_prev - the previous solved and stored G_n data
*/
	complex<double> *fx,*fy;
	fx = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	fy = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	init(fx,nx1*nx2);
	init(fy,nx1*nx2);
	mult(nx1*nx2,alpha1p,&fpn[nx1*nx2],fx);
	mult(nx1*nx2,one,fx,fy);

	mult(nx1*nx2,alpha2p,&fpn[nx1*nx2],fy);
	mult(nx1*nx2,one,fy,fy);

	for(int i=0;i<nx1;i++)
		for(int j=0;j<nx2;j++)
		{
			complex<double> *right1,*right2,*right3;
			right1 = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
			right2 = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
			right3 = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
			init(right1,3*nx2*nx1);
			init(right2,3*nx2*nx1);
			init(right3,3*nx2*nx1);

			int ind = i*nx2 + j;	
			int ind1 = nx1*nx2*((n-1)*nx1*nx2*9 + i*nx2*9 + j*9);
			int ind2 = nx1*nx2*(n*nx1*nx2*9 + i*nx2*9 + j*9);

			stress_n(&Un1[ind1],&Un1[ind2],alpha1p[ind],alpha2p[ind],betap1[ind],betap2[ind],lambda,mu,nx1,nx2,fx,fy,right1);
			stress_n(&Un2[ind1],&Un2[ind2],alpha1p[ind],alpha2p[ind],betap1[ind],betap2[ind],lambda,mu,nx1,nx2,fx,fy,right2);
			stress_n(&Un3[ind1],&Un3[ind2],alpha1p[ind],alpha2p[ind],betap1[ind],betap2[ind],lambda,mu,nx1,nx2,fx,fy,right3);

			complex<double> *temp111,*temp211,*temp311;
			temp111 = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
			temp211 = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
			temp311 = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
			init(temp111,3*nx2*nx1);
			init(temp211,3*nx2*nx1);
			init(temp311,3*nx2*nx1);
	




			delete right1,right2,right3;
			delete temp111,temp211,temp311;


		}

	delete fx,fy;

}



void solve_DNSO(int nx1,int nx2,int N,complex<double>* fpn,double* alpha1p,double* alpha2p,complex<double>* betap1,complex<double>* betap2,complex<double>* Un1,complex<double>* Un2,complex<double>* Un3,double lambda,double mu,complex<double>* Gn_m_mUm1,complex<double>* Gn_m_mUm2,complex<double>* Gn_m_mUm3)
{

/*
% solve_DNSO - Computes the Fourier coefficients
%     of the n-th Taylor term of the field in the upper domain
%     at the interface, U_n
%
% Inputs:
%
% . nx1 - Number of equally spaced gridpoints in x1
% . nx2 - Number of equally spaced gridpoints in x2
% . N - Number of perturbation orders
% . fpn - Fourier coefficients of f^n/n!
% . alpha1p - quasiwavenumbers: \alpha_1 + (2 \pi/d_1) p
% . alpha2p - quasiwavenumbers: \alpha_2 + (2 \pi/d_2) p
% . betap - \beta_p = \sqrt{ k^2 - \alpha_p^2 }
% . zetahat_n1 - Fourier coefficients of n-th Taylor term of
%    Dirichlet data
%
% Output:
%
% . Un - Fourier coefficients of U_n
% . Gn_m_mUm - Fourier coefficients of G_{n-m}[V_m]
%
% Note: Fourier coefficients of f^n/n! stored in fpn(:,n+1)
% Note: Fourier coefficients of U_n = (V_l)_n stored in Vln(:,n+1)
% Note: Fourier coefficients of G_{n-m}[U_m] stored in Gn_m_mUm(:,n-m+1,m+1)
*/
	complex<double> alpha1,alpha2,beta1,beta2,k1,k2;
	alpha1 = alpha1p[0];
	alpha2 = alpha2p[0];
	beta1 = betap1[0];
	beta2 = betap2[0];
	k1 = sqrt(pow(alpha1,2) + pow(alpha2,2) + pow(beta1,2));	
	k2 = sqrt(pow(alpha1,2) + pow(alpha2,2) + pow(beta2,2));	

/*
% Gn_m_mUm = G_{n-m}[U_m]
% Gn_m_mUm1.2.3 contains all orders of fourier expansion coefficients of G
% operator acting on all different wave numbers and different dimensions,
% here 3 would stands for different dimension in result, which is the
% result after the operation.
*/
	for(int n=0;n<N+1;n++)
	{
		if(n == 0)
			g0_engine(n,nx1,nx2,Un1,Un2,Un3,fpn,alpha1p,alpha2p,betap1,betap2,lambda,mu,Gn_m_mUm1,Gn_m_mUm2,Gn_m_mUm3); 

		else
			gn_engine(n,nx1,nx2,Un1,Un2,Un3,fpn,alpha1p,alpha2p,betap1,betap2,lambda,mu,Gn_m_mUm1,Gn_m_mUm2,Gn_m_mUm3); 





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
/*
	for(int i=0;i<nx1;i++)
		for(int l=0;l<nx2;l++)
		{
			xx1[nx2*i + l] = 0;
			xx2[nx2*i + l] = 0;
		}
*/
	init(xx1,nx1*nx2);
	init(xx2,nx1*nx2);

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
/*
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
*/
	init(pp1,nx1*nx2);
	init(pp2,nx1*nx2);
	init(alpha1p,nx1*nx2);
	init(alpha2p,nx1*nx2);
	init(betap1,nx1*nx2);
	init(betap2,nx1*nx2);

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
/*
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
*/

	init(f,nx1*nx2);
	init(fx1,nx1*nx2);
	init(fx1hat,nx1*nx2);
	init(fx2,nx1*nx2);
	init(fx2hat,nx1*nx2);
/*	
	for(int i=0;i<max(2,N+1);i++)
		for(int j=0;j<nx1;j++)
			for(int l=0;l<nx2;l++)
			{
				int k = i*nx1*nx2 + j*nx2 + l;
				fpn[k] = 0.0;
			}
*/
	init(fpn,max(2,N+1)*nx1*nx2);

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
	for(int i=0;i<nx1;i++)
		for(int j=0;j<nx2;j++)
		{
			int ind=i*9*nx2 + j*9;
			int ind1 = i*nx2 + j;
			N_inv[ind + 0] = alpha1p[ind1];	
			N_inv[ind + 1] = -betap2[ind1];	
			N_inv[ind + 2] = 0;	
			N_inv[ind + 3] = alpha2p[ind1];	
			N_inv[ind + 4] = 0;	
			N_inv[ind + 5] = -betap2[ind1];	
			N_inv[ind + 6] = betap1[ind1];	
			N_inv[ind + 7] = alpha1p[ind1];	
			N_inv[ind + 8] = alpha2p[ind1];	
			complex<double> *tmp[3];
			for(int i1 = 0;i1<3;i1++)
			{
				tmp[i1] = (complex<double>*)calloc(3,sizeof(complex<double>));
				for(int j1=0;j1<3;j1++)
				{
					int ind2 = i1*3 + j1;
					tmp[i1][j1] = N_inv[ind + ind2];
				}
			}
			
			complex<double> tmp1[3]; // store (1 0 0)
			complex<double> tmp2[3]; // store (0 1 0)
			complex<double> tmp3[3]; // store (0 0 1)
			init(tmp1,3); tmp1[0] = 1;
			init(tmp2,3); tmp2[1] = 1;
			init(tmp3,3); tmp3[2] = 1;
			int ind3 = i*3*nx2 + j*3;
			complexGaussElimination(tmp,tmp1,&a1[ind3],3);
			complexGaussElimination(tmp,tmp2,&a2[ind3],3);
			complexGaussElimination(tmp,tmp3,&a3[ind3],3);
			/*
			if(i==0&&j==0)
			{
				cout<<tmp1[0]<<" "<<tmp2[1]<<" "<<tmp3[2]<<endl;
				for(int i1=0;i1<3;i1++) 
					for(int j1=0;j1<3;j1++)
						cout<<tmp[i1][j1]<<" ";
				cout<<endl;
				for(int i1=0;i1<3;i1++)
					cout<<a1[ind3 + i1]<<" ";
				cout<<endl;
				
			}
			*/

		}

	// check point 3, finished checking: 2013/5/17
	
// Now let's set up dirichlet boundary condition
	complex<double> A,gbar;
	A = gbar = 0;

	complex<double> *Un1,*Un2,*Un3,*temp,*temp2;
	Un1 = (complex<double>*)calloc(9*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
	Un2 = (complex<double>*)calloc(9*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
	Un3 = (complex<double>*)calloc(9*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
//	gbar = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
//	A = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	temp = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	temp2 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	init(Un1,9*(N+1)*nx1*nx2*nx1*nx2);	
	init(Un2,9*(N+1)*nx1*nx2*nx1*nx2);	
	init(Un3,9*(N+1)*nx1*nx2*nx1*nx2);	
//	init(gbar,nx1*nx2);
//	init(A,nx1*nx2);
	init(temp,nx1*nx2);
	init(temp2,nx1*nx2);

	for(int n=0;n<=N;n++)
// (1 0 0) decompose to k^(1), f^(2), e^(2)
// first dimension of k^(1), f^(2) and e^(2)
// as to Un, first 3 stands for 3 types of wave p, s_1, s_2
// second 3 stands for 3 component of a wave: (1 0 0) (0 1 0) (0 0 1)
	{
		for(int m=0;m<nx1;m++)
			for(int k=0;k<nx2;k++)
			{
				complex<double> *base_hat;
				base_hat = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				init(base_hat,nx1*nx2);
				int ind = m*nx2 + k;
				base_hat[ind] = 1;
				complex<double> egbar = exp(one*betap1[ind]*gbar);
				complex<double> egbar1 = exp(one*betap2[ind]*gbar);
				// here we calculate temp and temp2, which indicates the 
           		// n-th order taylor  
           		// expansion of exp(1i*betap*g)  
				A = pow(one*betap1[ind],n)*egbar;
				mult(nx1*nx2,A,&fpn[n*nx1*nx2],temp);
/*				for(int m1=0;m1<nx1;m1++)
					for(int k1=0;k1<nx2;k1++)
					{
						int ind = m1*nx2 + k1;
						temp[ind] = A + fpn[n*nx1*nx2 + ind];
					}
*/
				A = pow(one*betap2[ind],n)*egbar1;
				mult(nx1*nx2,A,&fpn[n*nx1*nx2],temp2);
/*
				for(int m1=0;m1<nx1;m1++)
					for(int k1=0;k1<nx2;k1++)
					{
						int ind = m1*nx2 + k1;
						temp[ind] = A + fpn[n*nx1*nx2 + ind];
					}
*/
			
				for(int i=0;i<3;i++)
					for(int j=0;j<3;j++)
					{
//rule of indexing:
//(1)the slowest changing dimension is the perturbation order, eg dimension 3 in the matlab matrix, 
//(2) the last two dimension in the matrix should follow one by one then,eg dim 6 and 7.
//(3) then wave type, eg dim 4 
//(4) then the wave dimension, eg dim 5.
//(5) lastly the first two dimension of the matrix. 
						int ind = nx1*nx2*(n*9*nx1*nx2 + m*nx2*9 + k*9 + i*3 + j);
						int ind1 = 3*(m*nx2 + k) + i;
						int ind2 = 9*(m*nx2 + k) + j*3 + i;
						if(i == 0)
						{
							mult(nx1*nx2,a1[ind1]*N_inv[ind2],temp,&Un1[ind]);	
							newconv3d(nx1,nx2,base_hat,&Un1[ind],&Un1[ind]);	

							mult(nx1*nx2,a2[ind1]*N_inv[ind2],temp,&Un2[ind]);	
							newconv3d(nx1,nx2,base_hat,&Un2[ind],&Un2[ind]);	

							mult(nx1*nx2,a3[ind1]*N_inv[ind2],temp,&Un3[ind]);	
							newconv3d(nx1,nx2,base_hat,&Un3[ind],&Un3[ind]);	
						}
						else
						{
							mult(nx1*nx2,a1[ind1]*N_inv[ind2],temp2,&Un1[ind]);	
							newconv3d(nx1,nx2,base_hat,&Un1[ind],&Un1[ind]);	
							mult(nx1*nx2,a2[ind1]*N_inv[ind2],temp2,&Un2[ind]);	
							newconv3d(nx1,nx2,base_hat,&Un2[ind],&Un2[ind]);	
							mult(nx1*nx2,a3[ind1]*N_inv[ind2],temp2,&Un3[ind]);	
							newconv3d(nx1,nx2,base_hat,&Un3[ind],&Un3[ind]);	
						}

					}	
				delete base_hat;
			}

	}
// check point 4, finished checking

// Now here comes the main function
// First just do the 0-th order perturbation
	complex<double> *Gn_m_qUq1,*Gn_m_qUq2, *Gn_m_qUq3;
	Gn_m_qUq1 = (complex<double>*)calloc(3*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
	Gn_m_qUq2 = (complex<double>*)calloc(3*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
	Gn_m_qUq3 = (complex<double>*)calloc(3*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
	solve_DNSO(nx1,nx2,N,fpn,alpha1p,alpha2p,betap1,betap2,Un1,Un2,Un3,lambda,mu,Gn_m_qUq1,Gn_m_qUq2,Gn_m_qUq3);

	ofstream o1("Gn1.dat");
	int ind = nx1*nx2*(1*nx2*3 + 3*3 + 2);
	for(int i=0;i<nx1;i++)
	{
		for(int j=0;j<nx2;j++)
		{ 
			
			int k=i*nx2 + j;
			//k = 3*k;
			//for(int l=0;l<3;l++)
			o1<<setprecision(20)<<Gn_m_qUq3[k+ind]<<" ";
		}
		o1<<endl;
	}


	o1.close();	 
	delete pp1,pp2,alpha1p,alpha2p,betap1,betap2;
	delete fx1,fx2,fhat,fx2hat,fx1hat,fpn;
	delete N_inv, N_orig, a1,a2,a3;
	delete Un1, Un2, Un3, temp, temp2;
	delete Gn_m_qUq1,Gn_m_qUq2,Gn_m_qUq3;


	return 0;

}

