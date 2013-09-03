#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>




// "Minimal" ramdom number generator of Park and Miller with Bays-Durham shuffle and added safeguards 
// Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values)

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 3.0e-16
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */


double funk(double *v, int N){

  double **point, theta, phi, result, dif2;

  int i, j;


  point = (double **) malloc(2*N * sizeof(double *));

  for(i=0; i<(2*N); i++){

    point[i] = (double *) malloc(3 * sizeof(double));
  }


  for(i=0; i<N; i++){

    theta=v[2*i];
    phi=v[2*i+1];

    point[i][0] = sin(theta)*cos(phi);
    point[i][1] = sin(theta)*sin(phi);
    point[i][2] = cos(theta);

    point[N+i][0] = -point[i][0];
    point[N+i][1] = -point[i][1];
    point[N+i][2] = -point[i][2];

  }

/*   for(i=0; i<(2*N); i++){ */

/*     printf("point %i: %lf %lf %lf\n", i, point[i][0], point[i][1], point[i][2]); */
/*   } */



  result=0;

  for(i=0; i<(2*N); i++){
    for(j=0; j<(2*N); j++){

      if(j>i){

	dif2 = pow((point[i][0]-point[j][0]), 2.0)+pow((point[i][1]-point[j][1]), 2.0)+pow((point[i][2]-point[j][2]), 2.0);

	result += 1/dif2;
      }
    }
  }



  for(i=0; i<(2*N); i++){

    free(point[i]);
  }

  free(point);


  return result;

}



#define NRANSI


double amotry(double **p, double *y, double *psum, int ndim,
	      double (*funk)(double *, int), int ihi, double fac)
{
  int j, N;
  double fac1,fac2,ytry,*ptry;
 
  N=ndim/2;
 
  //ptry=vector(1,ndim);

  ptry = (double *)malloc(ndim * sizeof(double));
  
  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
  ytry=(*funk)(ptry, N);
  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j=0;j<ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  
  
  //	free_vector(ptry,1,ndim);
  
  free(ptry);
  
  return ytry;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */



void GET_PSUM(double **p, double *psum, int ndim, int mpts){
							
  int i, j;

  double sum;

  for (j=0;j<ndim;j++) {
		
    sum=0.0;
				
    for (i=0; i<mpts; i++) sum += p[i][j];			

    psum[j]=sum;
  
  }
}




#define NRANSI

#define NMAX 100000000

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}


void amoeba(double **p, double *y, int ndim, double ftol,
	    double (*funk)(double *, int), int *nfunk)
{

  double amotry(double **p, double *y, double *psum, int ndim,
		double (*funk)(double *, int), int ihi, double fac);

  int i,ihi,ilo,inhi,j,mpts=ndim+1, N=ndim/2;
  double rtol,sum,swap,ysave,ytry,*psum;
  
//  psum=vector(1,ndim);


  psum = (double *)malloc(ndim * sizeof(double));

  *nfunk=0;


  GET_PSUM(p, psum, ndim, mpts);

  
  for (;;) {
    ilo=0;
    ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
    for (i=0;i<mpts;i++) {
      if (y[i] <= y[ilo]) ilo=i;
      if (y[i] > y[ihi]) {
	inhi=ihi;
	ihi=i;
      } else if (y[i] > y[inhi] && i != ihi) inhi=i;
    }
    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
    if (rtol < ftol) {
      SWAP(y[0],y[ilo])
	for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i])
	  break;
      }
    if (*nfunk >= NMAX){
      
      printf("NMAX exceeded\n");
      exit(2);
    }
    
    *nfunk += 2;
    ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
    if (ytry <= y[ilo])
      ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
    else if (ytry >= y[inhi]) {
      ysave=y[ihi];
      ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
      if (ytry >= ysave) {
	for (i=0;i<mpts;i++) {
	  if (i != ilo) {
	    for (j=0;j<ndim;j++)
	      p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
	    y[i]=(*funk)(psum, N);
	  }
	}
	*nfunk += ndim;
	GET_PSUM(p, psum, ndim, mpts);
      }
    } else --(*nfunk);
  }
  // free_vector(psum,1,ndim);
  
  free(psum);
  
}
#undef SWAP
#undef GET_PSUM
#undef NMAX
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #!-"k$'%)]'V',4$. */




int main(int argc, char *argv[]){


  if(argc<3){

    printf("Usage: %s num_directions out_file\n", argv[0]);
    exit(1);
  }


  int N, ndim, i, j, *nfunk, index;

  double **p, *v, *y, theta, phi, ftol, min_y;

  FILE *output;

  long *idum;


  idum = (long *) malloc(sizeof(long));

  time_t ultime;
  time (&ultime); 

  idum[0]=-ultime;



  N=atoi(argv[1]);

  output = fopen(argv[2], "w");

  ndim = 2*N;

  p = (double **)malloc((ndim+1) * sizeof(double *));
  for(i=0; i<(ndim+1); i++){

    p[i] = (double *)malloc(ndim * sizeof(double));
  }

  v = (double *)malloc(ndim * sizeof(double));

  y = (double *)malloc((ndim+1) * sizeof(double));

  nfunk = (int *)malloc(1 * sizeof(int));

  for(i=0; i<(ndim+1); i++){

    for(j=0; j<N; j++){

      theta = ran1(idum)*M_PI;
      phi = ran1(idum)*2*M_PI;

 
      p[i][2*j] = theta;
      p[i][2*j+1] = phi;

    }
  }



  for(i=0; i<(ndim+1); i++){

    for(j=0; j<ndim; j++){

      v[j] = p[i][j]; 
    }

    y[i]=funk(v, N);

    //  printf("y[%i]=%lf\n", i, y[i]);

  }


  ftol=10e-10;



  amoeba(p, y, ndim, ftol, funk, nfunk);



  for(i=0; i<(ndim+1); i++){

    printf("arrangement %i:\n", i);

    for(j=0; j<N; j++){

      theta=p[i][2*j];
      phi=p[i][2*j+1];

      printf("direction %i: %lf, %lf %lf\n", j, sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    }
    printf("potential energy: %lf\n\n", y[i]);

  }

  /* Calculate the total gradient amplitude along each spatial coordinate for each arrangement of directions */
 
 double Gx, Gy, Gz, G_total[ndim+1];
 
 for(i=0; i<(ndim+1); i++){
    
    Gx=0.0;
    Gy=0.0;
    Gz=0.0;
    
    for(j=0; j<N; j++){
    
      theta=p[i][2*j];
      phi=p[i][2*j+1];
    
      Gx += fabs(sin(theta)*cos(phi));
      Gy += fabs(sin(theta)*sin(phi));
      Gz += fabs(cos(theta));
   }
   
   G_total[i]=Gx+Gy+Gz;
   
 }

 
 double G_min;

 G_min=G_total[0];
 index=0;

 for(i=1; i<(ndim+1); i++){
  
   if(G_total[i]<G_min){
   
     G_min=G_total[i];
     index=i;
   }

 }

 fprintf(output, "[directions=%i]\n", N);
 fprintf(output, "coordinatesystem[0]=prs\n");
 fprintf(output, "normalisation[0]=unity\n");

 for(j=0; j<N; j++){

    theta = p[index][2*j];
    phi = p[index][2*j+1];

    fprintf(output, "vector[%i]=(%lf, %lf, %lf)\n", j, sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    /*fprintf(output, "%lf, %lf, %lf\n", sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));*/

  }

  printf("directions written in %s\n", argv[2]);



  /*for(i=0; i<(ndim+1); i++){
  
    printf("G_total[%i]=%lf\n", i, G_total[i]);
    
  }

  printf("Minimum total gradient amplitude: index=%i, G_total=%lf\n", index, G_min);
*/

  free(nfunk);
  free(idum);
  free(v);
  free(y);

  for(i=0; i<(ndim+1); i++){

    free(p[i]);
  }

  free(p);


  return 0;
}
