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
	fftw_free(in);
	fftw_free(out);
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
	fftw_free(in);
	fftw_free(out);
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
//	f_temp = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
//	f_temp_hat = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
	f_temp = new complex<double>[nx1*nx2];
	f_temp_hat = new complex<double>[nx1*nx2];
	
//	f_temp = new complex<double>(nx1*nx2);
//	f_temp_hat = new complex<double>(nx1*nx2);
	init(f_temp,nx1*nx2);
	init(f_temp_hat,nx1*nx2);

	ofstream dd("test.out");
	for(int i=0;i<nx2*nx1;i++)
	{
		f_temp[i] = a*a1* exp(one*b_1*epsilon*f[i])
				  + b*b1* exp(one*b_2*epsilon*f[i]) 
				  + c*c1* exp(one*b_2*epsilon*f[i]);
		//cout<<"f_temp["<<i<<"]= "<<f_temp[i]<<endl;
		dd<<f_temp[i]<<" ";
/*
		
		cout<<"i = "<<i<<" "<< a*a1* exp(one*b_1*epsilon*f[i]);
		cout<<" a = "<<a<<" a1="<<a1<<" b_1 = "<<b_1;
		cout<<b*b1* exp(one*b_2*epsilon*f[i]);
		cout<<c*c1* exp(one*b_2*epsilon*f[i]);
		//if(radius(f_temp[i]) > 100) 
		cout<<"f_temp[i] = "<<f_temp[i]<<endl;
*/
	}
	
//	cout<<"Lover"<<endl;
/*
	cout<<"size of complex<double> is "<<sizeof(complex<double>)<<endl;
	cout<<"f_temp is "<<sizeof(*f_temp)<<" bytes!"<<endl;
	for(int i = 0;i<nx1;i++)
	{
		for(int j=0;j<nx2;j++)
		 	cout<<f_temp[i*nx2 + j]<<" ";
		cout<<endl;
	}
*/
	
	fftw_complex *in,*out;
	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(nx1*nx2));
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(nx1*nx2));

	c2fc(f_temp,in,nx1*nx2);

	fftw_plan p;
//	cout<<"Lover"<<endl;
	p = fftw_plan_dft_2d(nx1,nx2,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	
	fc2c(out,f_temp_hat,nx1*nx2);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	fft2(f_temp,nx1,nx2,f_temp_hat);
	newconv3d(nx1,nx2,f_temp_hat,base_hat,bc);
/*
*/
	delete f_temp,f_temp_hat;
	//free(f_temp);
	//free(f_temp_hat);
			
}

void direct_cal(int nx1,int nx2,complex<double> a1,complex<double> a2,complex<double> b1,complex<double> b2,complex<double> a,complex<double> b,complex<double> c,complex<double>* f,double epsilon,int m,int n,complex<double>* bc1,complex<double>* bc2,complex<double>* bc3)
{
	complex<double>* base_hat;
	base_hat = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
//	base_hat = new complex<double>(nx1*nx2);
	init(base_hat,nx1*nx2);
	base_hat[m*nx2 + n] = 1;
	direct_cal_deep(nx1,nx2,a,b,c,a1,-b2,0,b1,b2,f,epsilon,base_hat,bc1);

	delete base_hat;	
	//free(base_hat);
}

int main()
{

/*
% N - Number of Taylor orders
% nx1 - Number of equally spaced gridpoints on [0,d1]
% nx2 - Number of equally spaced gridpoints on [0,d2]
*/
	
	//ifstream o1("DNSO_cpp1.out");

	ifstream o1("DNSO_cpp.out");
	int N,nx1,nx2;
	double epsilon,d1,d2,dx1,dx2;
	o1>>N>>nx1>>nx2>>epsilon>>d1>>d2;
	cout<<"N = "<<N<<" nx1 = "<<nx1<<" nx2 = "<<nx2<<" epsilon = "<<epsilon<<" d1 = "<<d1<<" d2 = "<<d2<<endl;
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
	double *alpha1p,*alpha2p;///`j,*pp1,*pp2;
	complex<double> *betap1,*betap2;
	alpha1p = (double*)calloc(nx1*nx2,sizeof(double));
	alpha2p = (double*)calloc(nx1*nx2,sizeof(double));
/*
	pp1 = (double*)calloc(nx1*nx2,sizeof(double));
	pp2 = (double*)calloc(nx1*nx2,sizeof(double));
*/
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
//	for(int i=0;i<nx1*nx2;i++) cout<<"f["<<i<<"] = "<<f[i]<<endl;

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
	complex<double> tmp1[3]; // store (1 0 0)
	complex<double> tmp2[3]; // store (0 1 0)
	complex<double> tmp3[3]; // store (0 0 1)
	init(tmp1,3); tmp1[0] = 1;
	init(tmp2,3); tmp2[1] = 1;
	init(tmp3,3); tmp3[2] = 1;
	

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
//	cout<<"i = "<<i<<" j = "<<j<<" love"<<endl;
/*			
			complex<double> tmp1[3]; // store (1 0 0)
			complex<double> tmp2[3]; // store (0 1 0)
			complex<double> tmp3[3]; // store (0 0 1)
			init(tmp1,3); tmp1[0] = 1;
			init(tmp2,3); tmp2[1] = 1;
			init(tmp3,3); tmp3[2] = 1;
*/
			int ind3 = i*3*nx2 + j*3;
			complexGaussElimination(tmp,tmp1,&a1[ind3],3);
			complexGaussElimination(tmp,tmp2,&a2[ind3],3);
			complexGaussElimination(tmp,tmp3,&a3[ind3],3);
			
		//	for(int i1 = 0;i1<3;i1++)
		//		delete tmp[i1];
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

/*
	for(int i1=0;i1<1;i1++)
	{
		for(int j1=0;j1<2;j1++)
		{
			for(int k1 =0;k1<9;k1++)
			cout<<N_inv[9*(i1*nx2 + j1) + k1]<<" ";
			cout<<endl;
		}
		cout<<endl;
	}
*/

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
				//N_inv_test = new complex<double>*[3];
				ind = 9*(m*nx2 + n);

				for(int i1 = 0;i1<3;i1++)
				{
					N_inv_test[i1] = (complex<double>*)calloc(3,sizeof(complex<double>));
			//		N_inv_test[i1] = new complex<double>(3);
					for(int j1=0;j1<3;j1++)
					{
//	cout<<"i = "<<i1<<" j ="<<j1<<" love"<<endl;
						int ind2 = i1*3 + j1;
						N_inv_test[i1][j1] = N_inv[ind + ind2];
						//cout<<N_inv_test[i1][j1]<<" ";
					}
					//cout<<endl;
				}

				a1 = N_inv_test[0][0];
				a2 = N_inv_test[1][0];
				b1 = N_inv_test[2][0];
				b2 = -N_inv_test[0][1];

				complex<double> *bb,*cc;
				bb = (complex<double>*)calloc(3,sizeof(complex<double>));
				cc = (complex<double>*)calloc(3,sizeof(complex<double>));
				init(bb,3);
				init(cc,3);
				bb[dim] = 1;
			
				complexGaussElimination(N_inv_test,bb,cc,3);
				a = cc[0]; b = cc[1]; c = cc[2];


				complex<double> *bc1,*bc2,*bc3;
				bc1 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				bc2 = (complex<double>*)calloc(nx2*nx1,sizeof(complex<double>));
				bc3 = (complex<double>*)calloc(nx2*nx1,sizeof(complex<double>));
				
				direct_cal(nx1,nx2,a1,a2,b1,b2,a,b,c,f,epsilon,m,n,bc1,bc2,bc3);
				
				complex<double> *result,*base,*base_hat,*z1,*z2,*base_1,*base_2;
				result = (complex<double>*)calloc(nx1*nx2*3,sizeof(complex<double>));
				base_hat = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				base = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				base_1 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				base_2 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				z1 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				z2 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));

				init(result,nx1*nx2*3);
				init(base_hat,nx1*nx2);
				init(base,nx1*nx2);
				init(base_1,nx1*nx2);
				init(base_2,nx1*nx2);
				init(z2,nx1*nx2);
				init(z1,nx1*nx2);
				for(int order = 0;order <N+1;order++)	
					for(int k1=0;k1<nx1;k1++)
						for(int k2=0;k2<nx2;k2++)
						{
							int ind = nx1*nx2*3*(order*nx1*nx2 + k1*nx2 + k2);
							int ind1 = k1*nx2 + k2;
							for(int kk = 0;kk<nx1*nx2*3;kk++)
							result[kk] += ( bc1[ind1]*Gn_m_qUq1[ind + kk]
										+ bc2[ind1]*Gn_m_qUq2[ind + kk]
										+ bc3[ind1]*Gn_m_qUq3[ind + kk]
										)*pow(epsilon,order); 
						}
        			
				ifft2(result,nx1,nx2,result);
				ifft2(&result[nx1*nx2],nx1,nx2,&result[nx1*nx2]);
				ifft2(&result[2*nx1*nx2],nx1,nx2,&result[nx1*nx2*2]);

// here we need to check the real result	
				base_hat[m*nx2 + n] = 1;	
				ifft2(base_hat,nx1,nx2,base);
				for(int i1 = 0;i1<nx1*nx2;i1++)
				{
					z1[i1] = exp(one*b1*epsilon*f[i1]);
					z2[i1] = exp(one*b2*epsilon*f[i1]);
					base_1[i1] = z1[i1]*base[i1];
					base_2[i1] = z2[i1]*base[i1];
				} // in physical space
				
				complex<double> *p_wave,*s_wave,*wave_d1,*wave_d2,*wave_d3,
								*DIV,*u_partial_1_2,*u_partial_1_3,
								*u_partial_2_3,*u_3_partial_3,*G_result;

				p_wave = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
				s_wave = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));
				G_result = (complex<double>*)calloc(3*nx1*nx2,sizeof(complex<double>));

				wave_d1 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				wave_d2 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				wave_d3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				u_partial_1_2 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				u_partial_1_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				u_partial_2_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				u_3_partial_3 = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				DIV = (complex<double>*)calloc(nx1*nx2,sizeof(complex<double>));
				init(p_wave,3*nx1*nx2);
				init(s_wave,3*nx1*nx2);
				init(G_result,3*nx1*nx2);
				init(wave_d1,nx1*nx2);
				init(wave_d2,nx1*nx2);
				init(wave_d3,nx1*nx2);
				init(u_partial_1_2,nx1*nx2);
				init(u_partial_1_3,nx1*nx2);
				init(u_partial_2_3,nx1*nx2);
				init(u_3_partial_3,nx1*nx2);
				init(DIV,nx1*nx2);


				mult(nx1*nx2,a*a1,base_1,p_wave);
				mult(nx1*nx2,a*a2,base_1,&p_wave[nx1*nx2]);
				mult(nx1*nx2,a*b1,base_1,&p_wave[2*nx1*nx2]);

				mult(nx1*nx2,-b*b2,base_2,s_wave);
				mult(nx1*nx2,-c*b2,base_2,&s_wave[nx1*nx2]);
				mult(nx1*nx2,a1*b + c*a2,base_2,&s_wave[2*nx1*nx2]);

				for(int i1=0;i1<nx1*nx2;i1++)
				{
					wave_d1[i1] = p_wave[i1] + s_wave[i1];	
					wave_d2[i1] = p_wave[nx1*nx2 + i1] + s_wave[nx1*nx2 + i1];	
					wave_d3[i1] = p_wave[2*nx1*nx2+i1] + s_wave[2*nx1*nx2+i1];	
				}
// useful patterns, which are all in Phy-S
				mult(nx1*nx2,one*a*pow(k_1,2),base_1,DIV);
/*
				complex<double> aa = one*a*pow(k_1,2);
cout<<"love"<<endl;
					DIV[i1] = one*a*pow(k_1,2)*base_1[i1];
*/
				for(int i1 = 0;i1<nx1*nx2;i1++)
				{
					u_partial_1_2[i1] = one*(a1*wave_d2[i1] + a2*wave_d1[i1]);
					u_partial_1_3[i1] = one*(a1*wave_d3[i1] + 
										b1*p_wave[i1] +
										b2*s_wave[i1]);						
					u_partial_2_3[i1] = one*(a2*wave_d3[i1] +
										b1*p_wave[nx1*nx2 + i1] +
										b2*s_wave[nx1*nx2 + i1]);	
					u_3_partial_3[i1] = one*(b1*p_wave[2*nx1*nx2 + i1] +
										b2*s_wave[2*nx1*nx2 + i1]);
				}




				for(int i1=0;i1<3;i1++)
					delete N_inv_test[i1];
				delete N_inv_test;

				delete bb,cc;
				delete bc1,bc2,bc3;
				delete result,base,base_1,base_2,base_hat,z1,z2;
			}

		delete error;
	

	}

	delete alpha1p,alpha2p,betap1,betap2;
	delete fhat,f,fpn;
	delete N_inv, N_orig, a1,a2,a3;
	delete Un1, Un2, Un3;
	delete Gn_m_qUq1,Gn_m_qUq2,Gn_m_qUq3;
	return 0;

	




}


