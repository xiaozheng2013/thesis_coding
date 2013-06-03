/*
this is the c++ code for testing the acccuray of DNSO. 
*/
#include<iostream>
#include<fstream>
#include<iomanip>
#include<complex>
#include<fftw3.h>
#include<math.h>
#include<time.h>

#define pi 4.0*atan(1)
using namespace std;
complex<double> one(0,1);


void mult(int nx,complex<double> a, complex<double> *x, complex<double> *y)
{
	for(int i=0;i<nx;i++)
		y[i] = a*x[i];
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

double radius ( complex<double> c )
// Returns the radius of the complex number.
{
   double result,re,im;
   re = real(c);im=imag(c);
   result = re*re + im*im;
   return sqrt(result);
}


void mult(int nx,double* a, complex<double> *x, complex<double> *y)
{
	for(int i=0;i<nx;i++)
		y[i] = a[i]*x[i];
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


void direct_cal_deep(int nx1,int nx2,complex<double> a,complex<double> b,complex<double> c,complex<double> a1,complex<double> b1,complex<double> c1,complex<double> b_1,complex<double> b_2,complex<double>* f,double epsilon,complex<double>* base_hat,complex<double>* bc)
{
	complex<double> *f_temp,*f_temp_hat;
	f_temp = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	f_temp_hat = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));

	init(f_temp,nx1*nx2);
	init(f_temp_hat,nx1*nx2);

	for(int i=0;i<nx2*nx1;i++)
	{
		f_temp[i] = a*a1* exp(one*b_1*epsilon*f[i])
				  + b*b1* exp(one*b_2*epsilon*f[i]) 
				  + c*c1* exp(one*b_2*epsilon*f[i]);
	}
	fft2(f_temp,nx1,nx2,f_temp_hat);
	newconv3d(nx1,nx2,f_temp_hat,base_hat,bc);
	delete f_temp,f_temp_hat;
			
}

void direct_cal(int nx1,int nx2,complex<double> a1,complex<double> a2,complex<double> b1,complex<double> b2,complex<double> a,complex<double> b,complex<double> c,complex<double>* f,double epsilon,int m,int n,complex<double>* bc1,complex<double>* bc2,complex<double>* bc3)
{
	complex<double>* base_hat;
	base_hat = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	init(base_hat,nx1*nx2);
	base_hat[m*nx2 + n] = 1;
	direct_cal_deep(nx1,nx2,a,b,c,a1,-b2,0,b1,b2,f,epsilon,base_hat,bc1);
   

	delete base_hat;	
}

int main()
{

/*
% N - Number of Taylor orders
% nx1 - Number of equally spaced gridpoints on [0,d1]
% nx2 - Number of equally spaced gridpoints on [0,d2]
*/
	
	ifstream o1("DNSO_cpp1.out");
	int N,nx1,nx2;
	double epsilon,d1,d2,dx1,dx2;
	o1>>N>>nx1>>nx2>>epsilon>>d1>>d2;
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
	for(int i=0;i<nx1*nx2;i++)	o1>>alpha1p[i];
	for(int i=0;i<nx1*nx2;i++)	o1>>alpha2p[i];
	for(int i=0;i<nx1*nx2;i++)	o1>>betap1[i];
	for(int i=0;i<nx1*nx2;i++)	o1>>betap2[i];

// set up the powers of f
	
	complex<double> *fpn,*f,*fhat;

	fhat = (complex<double> *)calloc(nx1*nx2,sizeof(complex<double>));
	f = (complex<double> *)calloc(nx1*nx2,sizeof(complex<double>));

	init(f,nx1*nx2);
	init(fhat,nx1*nx2);
	fhatsetup3d(nx1,nx2,fhat);
	ifft2(fhat,nx1,nx2,f);

	fpn = (complex<double> *)calloc(nx1*nx2*max(2,N+1),sizeof(complex<double>));
	for(int i=0;i<nx2*nx1*max(2,N+1);i++) o1>>fpn[i];

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

		complex<double> *Un1,*Un2,*Un3;
		Un1 = (complex<double>*)calloc(9*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
		Un2 = (complex<double>*)calloc(9*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
		Un3 = (complex<double>*)calloc(9*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
		for(int i=0;i<(N+1)*9*nx1*nx2*nx1*nx2;i++)	o1>>Un1[i];
		for(int i=0;i<(N+1)*9*nx1*nx2*nx1*nx2;i++)	o1>>Un2[i];
		for(int i=0;i<(N+1)*9*nx1*nx2*nx1*nx2;i++)	o1>>Un3[i];

		complex<double> *Gn_m_qUq1,*Gn_m_qUq2, *Gn_m_qUq3;
		Gn_m_qUq1 = (complex<double>*)calloc(3*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
		Gn_m_qUq2 = (complex<double>*)calloc(3*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
		Gn_m_qUq3 = (complex<double>*)calloc(3*(N+1)*nx1*nx2*nx1*nx2,sizeof(complex<double>));
		for(int i=0;i<(N+1)*3*nx1*nx2*nx1*nx2;i++)	o1>>Gn_m_qUq1[i];
		for(int i=0;i<(N+1)*3*nx1*nx2*nx1*nx2;i++)	o1>>Gn_m_qUq2[i];
		for(int i=0;i<(N+1)*3*nx1*nx2*nx1*nx2;i++)	o1>>Gn_m_qUq3[i];
	
	for(int dim = 0;dim<3;dim++)
	{
		complex<double> *error;
		error = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
		complex<double> a1,a2,b1,b2,a,b,c;
		int ind;
		for(int m=0;m<nx1;m++)
			for(int n =0;n<nx2;n++)
			{
				complex<double> **N_inv_test;
				N_inv_test =(complex<double>**)calloc(3,sizeof(complex<double>*));
				ind = nx1*nx2*(m*nx2 + n);

				for(int i1 = 0;i1<3;i1++)
				{
					N_inv_test[i1] = (complex<double>*)calloc(3,sizeof(complex<double>));
					for(int j1=0;j1<3;j1++)
					{
						int ind2 = i1*3 + j1;
						N_inv_test[i1][j1] = N_inv[ind + ind2];
					}
				}

						

				a1 = N_inv_test[0][0];
				a2 = N_inv_test[1][0];
				b1 = N_inv_test[2][0];
				b2 = -N_inv_test[0][1];

				complex<double> *bb,*cc;
				bb = (complex<double>*)calloc(3,sizeof(complex<double>));
				cc = (complex<double>*)calloc(3,sizeof(complex<double>));
				init(bb,9);
				init(cc,9);
				bb[dim] = 1;
			
				complexGaussElimination(N_inv_test,bb,cc,3);
				a = cc[0]; b = cc[1]; c = cc[2];

				complex<double> *bc1,*bc2,*bc3;
				bc1 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				bc2 = (complex<double>*)calloc(nx2*nx1,sizeof(complex<double>));
				bc3 = (complex<double>*)calloc(nx2*nx1,sizeof(complex<double>));
				
				//direct_cal(nx1,nx2,a1,a2,b1,b2,a,b,c,f,epsilon,m,n,bc1,bc2,bc3);
        				
	/*	
				for(int i1=0;i1<3;i1++)
					delete N_inv_test[i1];
				delete[] N_inv_test;
*/
				delete bb,cc;
				delete bc1,bc2,bc3;

			}

		delete error;
	

	}
	




}


