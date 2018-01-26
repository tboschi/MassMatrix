#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>

typedef struct OSC_PARAMS { double theta12; double theta13; double theta23; double delta; double alpha21; double alpha31; double Dm21; double Dm31; double m3; double m2; double m1; } OSC_PARAMS;

int sinsq(double complex U[], double * a)
{
	double sinsqth13 = pow(cabs(U[2]),2.0);
	if(sinsqth13==1.0)
	{// THIS EXCEPTION NEEDS TO BE DEALT WITH.

		//printf("BAD ANGLE!\n"); 
		//return 1;
		//My logic is, although finding the angles is hard in this case. They are all ruled out as theta13=90. 

		a[1] = 100.0;
		a[0] = 100.0;//pow(cabs(U[3]),2.0)+pow(cabs(U[6]),2.0);
		a[2] = 100.0;
		a[3] = 100.0;
		a[4] = 100.0;
		a[5] = 100.0;
	}
	else
	{
		a[1] = sinsqth13;
		a[0] = pow(cabs(U[1]),2.0)/(1.0-sinsqth13);
		a[2] = pow(cabs(U[5]),2.0)/(1.0-sinsqth13);
		a[3] = carg(  -U[0]*conj(U[2]*U[3])*U[5] - (1-a[0])*(1-a[1])*a[1]*a[2]  );
		a[4] = (carg(U[1]*U[1]*conj(U[0]*U[0]))); // \alpha_{21}

		a[5] = carg(U[2]*U[2]*conj(U[0]*U[0])) + 2.0*a[3];

		if(a[5]>M_PI+1e-3){ a[5]-=2.0*M_PI;} 
		else if(a[5]<-M_PI-1e-3){ a[5]+=2.0*M_PI;}

		if(fabs(a[5]-M_PI)<1e-3){a[5]=-M_PI;}
 
//		printf("### a[5] = %.5lf\n", a[5]);

	}

return 0;
}

double find_PMNS(double complex M[], OSC_PARAMS * output)
{
  

  double complex PMNS[] = { 1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 
  			    0.0 + 0.0*I, 2.0 + 0.0*I, 9.0 + 0.0*I, 
			    0.0 + 0.0*I, 5.0 + 0.0*I, 3.0 + 0.0*I };

  double complex TEMP[] = { 1.0 + 0.0*I, 10.0 + 0.0*I, 3.0 + 0.0*I, 
  			    0.0 + 0.0*I, 2.0 + 0.0*I, 100.0 + 0.0*I, 
			    1.0 + 0.0*I, 5.0 + 0.0*I, 3.0 + 0.0*I };

  
  double complex beta = 0.0 + 0.0*I; 
  double complex alpha = 1.0 + 0.0*I;

  cblas_zgemm (CblasRowMajor, CblasConjTrans, CblasNoTrans, 3, 3, 3, &alpha, M, 3, M, 3, &beta, TEMP, 3);

  gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc(3);
  gsl_vector *eval = gsl_vector_alloc(3);
  gsl_matrix_complex  *evec = gsl_matrix_complex_alloc(3, 3);
  
  gsl_matrix_complex *mat = gsl_matrix_complex_alloc(3, 3);

  int x1,y1;
  for(x1=0;x1<3;x1+=1)
  {
   for(y1=0;y1<3;y1+=1)
   { 
 	gsl_matrix_complex_set(mat, x1, y1, gsl_complex_rect(creal(TEMP[x1*3+y1]),cimag(TEMP[x1*3+y1])));
   }
  }

  gsl_eigen_hermv(mat, eval, evec, w);
  gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

  int i,j;
  int n = 0;
  for(i=0;i<3;i++)
  {
		for(j=0;j<3;j++)
		{
			gsl_vector_complex_view outcol = gsl_matrix_complex_column (evec, i);
			gsl_complex z = gsl_vector_complex_get (&outcol.vector, j);
			PMNS[j*3+n] = GSL_REAL(z) + GSL_IMAG(z)*I;
			printf(" %g + %gi, ",  GSL_REAL(z), GSL_IMAG(z));
		}
		n+=1;
		printf("\n");	
  }

//printf("%.3lf %.3lf %.3lf\n\n", gsl_vector_get(eval,0), gsl_vector_get(eval,1), gsl_vector_get(eval,2)); 

   output->Dm21 = gsl_vector_get(eval,1) -  gsl_vector_get(eval,0); 
   output->Dm31 = gsl_vector_get(eval,2) -  gsl_vector_get(eval,0); 
   output->m3 = sqrt(gsl_vector_get(eval,2));
   output->m2 = sqrt(gsl_vector_get(eval,1));
   output->m1 = sqrt(gsl_vector_get(eval,0));

  double angles[6];

  sinsq(PMNS,angles);

  double theta12, theta13, theta23;
  if(angles[0]==100.0) { theta12 = 1000.0; }
  else {  theta12 = (180.0/M_PI)*asin(sqrt(angles[0])); }
  if(angles[1]==100.0) { theta13 = 1000.0; }
  else { theta13 = (180.0/M_PI)*asin(sqrt(angles[1])); }
  if(angles[2]==100.0) { theta23 = 1000.0; }
  else { theta23 = (180.0/M_PI)*asin(sqrt(angles[2])); }
	
  if(fabs(angles[0]-1.0)<1e-7){theta12 = 90.0;}
  if(fabs(angles[1]-1.0)<1e-7){theta13 = 90.0;}
  if(fabs(angles[2]-1.0)<1e-7){theta23 = 90.0;}

//  printf("%.3lf %.3lf %.3lf\n", theta12, theta13, theta23);

  output->theta12 = theta12;
  output->theta13 = theta13;
  output->theta23 = theta23;
  output->delta   = angles[3];
  output->alpha21 = angles[4];
  output->alpha31 = angles[5];

//  printf("%.3lf %.3lf %.3lf\n", output->theta12, output->theta13, output->theta23);

return 0;  
}

double find_PMNS_8dim(double complex M[], OSC_PARAMS * output)
{
  
  double complex PMNS[] = { 
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I};

  double complex TEMP[] = { 
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I,
	1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I,  0.0 + 0.0*I};

  double complex beta = 0.0 + 0.0*I; 
  double complex alpha = 1.0 + 0.0*I;

  cblas_zgemm (CblasRowMajor, CblasConjTrans, CblasNoTrans, 8, 8, 8, &alpha, M, 8, M, 8, &beta, TEMP, 8);

  gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc(8);
  gsl_vector *eval = gsl_vector_alloc(8);
  gsl_matrix_complex  *evec = gsl_matrix_complex_alloc(8, 8);
  
  gsl_matrix_complex *mat = gsl_matrix_complex_alloc(8, 8);

  int x1,y1;
  for(x1=0;x1<8;x1+=1)
  {
   for(y1=0;y1<8;y1+=1)
   { 
 	gsl_matrix_complex_set(mat, x1, y1, gsl_complex_rect(creal(TEMP[x1*8+y1]),cimag(TEMP[x1*8+y1])));
   }
  }

  gsl_eigen_hermv(mat, eval, evec, w);
  gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  int i,j;
  int n = 0;

  for(i=0;i<8;i++)
  {
		for(j=0;j<8;j++)
		{
			gsl_vector_complex_view outcol = gsl_matrix_complex_column (evec, i);
			gsl_complex z = gsl_vector_complex_get (&outcol.vector, j);
			PMNS[j*8+n] = GSL_REAL(z) + GSL_IMAG(z)*I;
			//printf(" %g + %gi, ",  GSL_REAL(z), GSL_IMAG(z));
			printf(" %g , ",  GSL_REAL(z));
		}
		n+=1;
		printf("\n");	
  }

  printf("\nMasses:\n");	
  printf("Dm21 = %g eV^2\n",  fabs(gsl_vector_get(eval, 1))-fabs(gsl_vector_get(eval,0)));
  printf("Dm31 = %g eV^2\n",  fabs(gsl_vector_get(eval, 2))-fabs(gsl_vector_get(eval,0)));
  printf("m_1 = %g eV\n",  sqrt(fabs(gsl_vector_get(eval, 0))));
  printf("m_2 = %g eV\n",  sqrt(fabs(gsl_vector_get(eval, 1))));
  printf("m_3 = %g eV\n",  sqrt(fabs(gsl_vector_get(eval, 2))));
  printf("m_s = %g eV\n",  sqrt(fabs(gsl_vector_get(eval, 3))));
  printf("m_H1 = %g eV\n",  sqrt(fabs(gsl_vector_get(eval, 4))));
  printf("m_H2 = %g eV\n",  sqrt(fabs(gsl_vector_get(eval, 5))));
  printf("m_H3 = %g eV\n",  sqrt(fabs(gsl_vector_get(eval, 6))));
  printf("m_H4 = %g eV\n",  sqrt(fabs(gsl_vector_get(eval, 7))));

//
//   output->Dm21 = gsl_vector_get(eval,1) -  gsl_vector_get(eval,0); 
//   output->Dm31 = gsl_vector_get(eval,2) -  gsl_vector_get(eval,0); 
//   output->m3 = sqrt(gsl_vector_get(eval,2));
//   output->m2 = sqrt(gsl_vector_get(eval,1));
//   output->m1 = sqrt(gsl_vector_get(eval,0));
//
//  double angles[6];
//
//  sinsq(PMNS,angles);
//
//  double theta12, theta13, theta23;
//  if(angles[0]==100.0) { theta12 = 1000.0; }
//  else {  theta12 = (180.0/M_PI)*asin(sqrt(angles[0])); }
//  if(angles[1]==100.0) { theta13 = 1000.0; }
//  else { theta13 = (180.0/M_PI)*asin(sqrt(angles[1])); }
//  if(angles[2]==100.0) { theta23 = 1000.0; }
//  else { theta23 = (180.0/M_PI)*asin(sqrt(angles[2])); }
//	
//  if(fabs(angles[0]-1.0)<1e-7){theta12 = 90.0;}
//  if(fabs(angles[1]-1.0)<1e-7){theta13 = 90.0;}
//  if(fabs(angles[2]-1.0)<1e-7){theta23 = 90.0;}
//
////  printf("%.3lf %.3lf %.3lf\n", theta12, theta13, theta23);
//
//  output->theta12 = theta12;
//  output->theta13 = theta13;
//  output->theta23 = theta23;
//  output->delta   = angles[3];
//  output->alpha21 = angles[4];
//  output->alpha31 = angles[5];
//
////  printf("%.3lf %.3lf %.3lf\n", output->theta12, output->theta13, output->theta23);
//
return 0;  
}



int main(int argc, char* argv[])
{ 

  static OSC_PARAMS out; 

//  double complex M[] = { 1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 
//  			 0.0 + 0.0*I, 3.0 + 0.0*I, 0.0 + 0.0*I, 
//			 0.0 + 0.0*I, 0.0 + 0.0*I, 3.0 + 0.0*I };
//
//  find_PMNS(M,&out);
 
  double complex m11 = 1.0e6 + 0.0*I;
  double complex m12 = 2.0e6 + 0.0*I;
  double complex m21 = 3.0e6 + 0.0*I;
  double complex m22 = 4.0e6 + 0.0*I;
  double complex m31 = 5.0e6 + 0.0*I;
  double complex m32 = 6.0e6 + 0.0*I;

  double complex L11 = 1e9 + 0.0*I;
  double complex L12 = 3e9 + 0.0*I;
  double complex L13 = 5e9 + 0.0*I;
  double complex L21 = 6e9 + 0.0*I;
  double complex L22 = 2e9 + 0.0*I;
  double complex L23 = 8e9 + 0.0*I;

  double complex mu11 = 4e2 + 0.0*I;
  double complex mu12 = 3e2 + 0.0*I;
  double complex mu13 = 5e2 + 0.0*I;
  double complex mu22 = 2e2 + 0.0*I;
  double complex mu23 = 6e2 + 0.0*I;
  double complex mu33 = 1e2 + 0.0*I;

  double complex M8[] = { 
	0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, m11,         m12,         0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 
	0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, m21,         m22,         0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 
	0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, m31,         m32,         0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 
	m11,         m21,         m31,         0.0 + 0.0*I, 0.0 + 0.0*I, L11,         L12,         L13,  
	m12,         m22,         m32,         0.0 + 0.0*I, 0.0 + 0.0*I, L21,         L22,         L23, 
	0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, L11,         L21,         mu11,        mu12,        mu13, 
	0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, L12,         L22,         mu12,        mu22,        mu23, 
	0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, L13,         L23,         mu13,        mu23,        mu33}; 

	int i = 0, j = 0;
  for(i=0;i<8;i++)
  {
		for(j=0;j<8;j++)
		{
			printf(" %g + i%g, ", creal(M8[8*i+j]), cimag(M8[8*i+j]));
		}
		printf("\n");	
  }

find_PMNS_8dim(M8,&out);
 
return 0;
}

