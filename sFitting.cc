/* sFitting.cc
   Photon Spectrum Fitting for Verify
10 MAY 2006 JOK created
*/

/* *********************************************************************** */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>             // for INT_MAX
#include <math.h>               // for pow
#include <string.h>             // for strcpy
#include <malloc.h>             // For memory allocation functions
#include <ctype.h>

using namespace std;
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) > (b) ? (b) : (a))


#define PASS 				0
#define FAIL				-1
#define MAX_STR_LEN		256

#define NRANSI
#include "nr.h"
#include "nrutil.h"

// boolean type
enum boolean {FALSE, TRUE};

// constants
const double MINUS_ONE = -1.0;
const double ZERO      =  0.0;
const double HALF      =  0.5;
const double ONE       =  1.0;
const double TWO       =  2.0;
const double SQRT2     = sqrt(2.0);
const double TWOSQRT2  = 2.0*SQRT2;

const double PI       = 3.14159265358979324;   // PI
const double TWO_PI   = 6.28318530717958648;   // 2*PI


void usage(){
  char *FuncName="sFitting";
  printf("\n Usage: %s -i mData -o oData -l Lval -b Bval -Emax Emax [-Emin Emin] [-norm norm] [-d] [-h] \n",FuncName);
  printf("\n        mData is a two columed file consisting of energy and prob.");
  printf("\n");
}

double P(double E, double L, double B, double Emax, double Emin, double N) {
	return N*(1-exp(-L*E))*(exp(-B*E) - exp(-B*Emax));
}

double dPdL(double E, double L, double B, double Emax, double Emin, double N) {
	return N*(0.0+E*exp(-L*E))*(exp(-B*E) - exp(-B*Emax));
}

double dPdB(double E, double L, double B, double Emax, double Emin, double N) {
	return N*(1-exp(-L*E))*(-E*exp(-B*E) + Emax* exp(-B*Emax));
}

double dPdEmax(double E, double L, double B, double Emax, double Emin, double N) {
	return N*(1-exp(-L*E))* B*exp(-B*Emax);
}

double dPdN(double E, double L, double B, double Emax, double Emin, double N) {
	return (1-exp(-L*E))*(exp(-B*E) - exp(-B*Emax));
}

// calculate dose and derivatives
void SpectCalc(float x, float a[], float *y, float dyda[], int na)
{
   // get the parameters
   if (na != 4) {
		 printf ("ERROR: wrong number of fit parameters in dcal\n");
		 exit (FAIL);
	}
	
   double L    = a[1];
	double B    = a[2];
	double Emax = a[3];
	double N    = a[4];
	
   double z = x;
	
	*y      = P(z, L, B, Emax, 0.25, N);
        dyda[1] = dPdL(z, L, B, Emax, 0.25, N);
	dyda[2] = dPdB(z, L, B, Emax, 0.25, N);
	dyda[3] = dPdEmax(z, L, B, Emax, 0.25, N);	
	dyda[4] = dPdN(z, L, B, Emax, 0.25, N);	
}


/* *********************************************************************** */
int main(int argc, char *argv[]){


	double Emax = 6.0;   int fitEmax = 1;
	double Emin = 0.25;  int fitEmin = 1;
	double L = 2.5;      int fitL = 1;
	double B = 0.45;     int fitB = 1;
	double N = 1500.0;   int fitN = 1;

  	/* Parse Command Line Arguments */
	char iFileName[MAX_STR_LEN];    // Input File Name
	char oFileName[MAX_STR_LEN];    // Output File Name

 	// Variable declaration for Flags
  	int   iFlag=0, oFlag=0 ;    // Flags for input and output files
	int   dFlag=0 ;             // Flags for debug
	int   xFlag=0, yFlag=0, zFlag=0;  // Flags for Coll X and Y and SAD Setting

  	// Get arguments from command line ----------------------------------------
  	iFileName[0] = 0;  // Clear the filename
	oFileName[0] = 0;  // Clear the filename
	
	float iEmax = (float) Emax;
	float iEmin = (float) Emin;
	float iL = (float) L;
	float iB = (float) B;	
	float iN = (float) N;

		
	
	// Check Arguments	
  	for(int iArg=0; iArg < argc; iArg++){
    	if(iArg < argc-1){
        	if( strcmp(argv[iArg],"-i") ==0 || strcmp(argv[iArg],"-input") == 0 ){
	    		iArg++;
	    		strcpy(iFileName,argv[iArg]);
				iFlag = 1;
	  		}
        	if( strcmp(argv[iArg],"-o") ==0 || strcmp(argv[iArg],"-output") == 0 ){
	    		iArg++;
	    		strcpy(oFileName,argv[iArg]);
				oFlag = 1;
	  		}
        	if( strcmp(argv[iArg],"-L") ==0 || strcmp(argv[iArg],"-l") == 0 ){
	    		iArg++;
	    		sscanf(argv[iArg], "%f", &iL);
				L = (double) iL;
				fitL = 0;
	  		}
        	if( strcmp(argv[iArg],"-B") ==0 || strcmp(argv[iArg],"-b") == 0 ){
	    		iArg++;
	    		sscanf(argv[iArg], "%f", &iB);
				B = (double) iB;
				fitB = 0;
	  		}
        	if( strcmp(argv[iArg],"-Emin") ==0 || strcmp(argv[iArg],"-emin") == 0 ){
	    		iArg++;
	    		sscanf(argv[iArg], "%f", &iEmin);
				Emin = (double) iEmin;
				fitEmin = 0;
	  		}
        	if( strcmp(argv[iArg],"-Emax") ==0 || strcmp(argv[iArg],"-emax") == 0 ){
	    		iArg++;
	    		sscanf(argv[iArg], "%f", &iEmax);
				Emax = (double) iEmax;
				fitEmax = 0;
	  		}
        	if( strcmp(argv[iArg],"-N") ==0 || strcmp(argv[iArg],"-norm") == 0 ){
	    		iArg++;
	    		sscanf(argv[iArg], "%f", &iN);
				N = (double) iN;
				fitN = 0;
	  		}
        	if( strcmp(argv[iArg],"-d") ==0 || strcmp(argv[iArg],"-debug") == 0 ){
	    		iArg++;
	    		sscanf(argv[iArg], "%d", &dFlag);
				printf ("# Debug Flag = %d\n",dFlag);
	  		}
      }
    	if(strcmp("-h", argv[iArg]) == 0 || strcmp("-help", argv[iArg]) == 0 ) {
        	usage();
        	return(PASS);
      }
  	}

   if(dFlag != 0) { 
      if (iFlag == 1) printf ("# Input File = %s\n", iFileName);
      if (oFlag == 1) printf ("# Output File = %s\n", oFileName);
   }
	
   /* Read Spectrum File -------------------------------------------------------- */	
 
	if(iFileName[0] == 0) {
      printf("\n ERROR: Must give INPUT filename\n");
      usage();
      return(FAIL);
	}
/* ===	
	if(oFileName[0] == 0) {
      printf("\n ERROR: Must give OUTPUT filename\n");
      usage();
      return(FAIL);
	}
=== */

	FILE *istrm;
	if( (istrm = fopen(iFileName,"r") ) == NULL ) {
		printf("\n ERROR: opening file %s \n", iFileName);
		return(FAIL);
	}
	rewind(istrm);

	/* === Read in the spectrum file === */
	int nData = 0; // number of data points 
	
	// fscanf (istrm, "# %d", &nData);
	// if (dFlag == 1) printf ("#  nData = %d\n", nData);
	int iLines = 0;
	char sLine[256];
	while (!feof(istrm)) {
		sLine[0] = 0;
		fgets(sLine, sizeof(sLine), istrm);
		// printf ("%d\n", strlen(sLine));
		if (strlen(sLine) > 1) iLines++;
	}
	printf ("iLines = %d\n", iLines);
	nData = iLines;
	rewind(istrm);
	
	/* === Memory Allocation for Dair === */  
	double *z = (double *) calloc (nData, sizeof(double));  // Energy
	double *Dm = (double *) calloc (nData, sizeof(double)); // Measured Spectrum
	double *Dc = (double *) calloc (nData, sizeof(double)); // Calculated Spectrum
		
	if (dFlag == 1) printf ("# Start reading measured spectrum --- \n");

	double peakEm = 0.0;
	double peakPm = 0.0;
	double zMax = 0.0;		
	for (int i=0;i<nData;i++) {
		float zValue, dValue;
		fscanf (istrm, "%e %e", &zValue, &dValue);
		z[i] = (double)zValue;
		Dm[i] = (double)dValue;

		if(Dm[i] >= peakPm) {
			 peakPm = Dm[i];
			 peakEm = z[i];
		}
		if(z[i] >= zMax) zMax = z[i];
	}	

	double peakEc = 0.0;
	double peakPc = 0.0;	
	// Emax =zMax;
	// N = 1.0;
	for (int i=0;i<nData;i++) {
		Dc[i] = P(z[i], L, B, Emax, Emin, N);
		// printf ("Dc[i]= %e\n", Dc[i]);
		if(Dc[i] >= peakPc) {
			 peakPc = Dc[i];
			 peakEc = z[i];
		}
	}	
		
	for (int i=0;i<nData;i++) {
		Dm[i] = Dm[i]/peakPm;
		Dc[i] = Dc[i]/peakPc;
		if (dFlag == 1) printf ("%f %f %f\n", z[i], Dm[i], Dc[i]);	
	}	

   const int  MZ = 2000;        // maximum number of z points

   double     dmeas[MZ];        // measured spectrum probability
   double     zmeas[MZ];        // measured sepctrum energy
   int        nzmeas;           // number of z points for measured curve

  // parameters and memory for the routines from Numerical Recipes (NR)
   const long NPT = MZ;            // maximum number of measured data points
   const long MA  = 4;             // number of fit parameters
   float      *x,*y,*sig;          // data points and variance
   int        *ia;
   float      *a;
   float      **covar,**alpha;
   float      chisq,chisq_old,alamda;

   // this is also NR stuff
   x   = vector(1,NPT);
   y   = vector(1,NPT);
   sig = vector(1,NPT);
   ia  = ivector(1,MA);
   a   = vector(1,MA);
   covar=matrix(1,MA,1,MA);
   alpha=matrix(1,MA,1,MA);

   // number of data points (size of x,y,sig)
   int ndata = 0;

   // normalize sigma
   double onesig2 = ZERO;
	
   // read z values and spectrum data
   // get spectrum
   int    k = 0;
   double dtmp;
   double ztmp = -1.0e10;
   double zold = ztmp;
   double dmeas_max = 0.0;
	boolean next = TRUE;
	
   next = TRUE;
   while (next) {
      zold = ztmp;
		ztmp = z[k];
		dtmp = Dm[k];
      if ((ztmp > zold) && (k<MZ)){
         zmeas[k] = ztmp;
         dmeas[k] = dtmp;
         if (dtmp > dmeas_max) dmeas_max = dtmp;

         double sigma = ONE;
         if (k>0) {
            double delx = fabs(ztmp-zmeas[k-1]);
            if (delx > 0.001) sigma = fabs((dtmp-dmeas[k-1])/delx);
            if (sigma < ONE)  sigma = ONE;
         }
         ++k;

         ++ndata;
         x[ndata] = ztmp;
         y[ndata] = dtmp;

         sig[ndata] = sigma; if (ndata == 2) sig[1] = sigma;
         if (ndata <= 2) onesig2  = TWO/(sigma*sigma);
         else            onesig2 += ONE/(sigma*sigma);
      } else {
         next = FALSE;
      }
   }
   nzmeas = k;

   // normalize measured data to 1
   double scale = 1.0/dmeas_max;
   for (k=0; k<nzmeas; ++k) dmeas[k] *= scale;
   dmeas_max = 1.0;
		
   // normalize standard deviation
   float onesig = sqrt(onesig2);
   for (int i=1; i<=ndata; ++i){
     sig[i] *= onesig;
   }	
	
	a[1] = L;     ia[1] = fitL;
	a[2] = B;     ia[2] = fitB;
	a[3] = Emax;  ia[3] = fitEmax;
	a[4] = N;     ia[4] = fitN;
	
   alamda = MINUS_ONE;
   mrqmin(x,y,sig,ndata,a,ia,MA,covar,alpha,&chisq,SpectCalc,&alamda);
	
   // perform the iterations
   int niter = 1;
   int itest = 0;
   do {
		printf ("#  Iteration: %d  chi: %f  alamda: %f \n", niter, sqrt(chisq), alamda);
		printf ("#  lval: %f  bval: %f  Emin: %f  Emax: %f  Norm: %f\n\n",
				     a[1], a[2], Emin, a[3], a[4]);

      ++niter;
      chisq_old = chisq;
      mrqmin(x,y,sig,ndata,a,ia,MA,covar,alpha,&chisq,SpectCalc,&alamda);
      if (chisq > chisq_old) {
         itest = 0;
      }
      else if (fabs(chisq_old-chisq) < 0.001) {
         ++itest;
      }
   } while (itest < 100 && alamda < 10000);
	
		
	printf ("# initial L  = %f ---> final L  = %f\n", L,  a[1]);
	printf ("# initial B = %f ---> final B = %f\n", B, a[2]);
	printf ("# initial Emax = %f ---> final Emax = %f\n", Emax, a[3]);
	printf ("# initial N = %f ---> final N = %f\n", N, a[4]);	

	if (dFlag == 1) printf ("# Start reading measured PDD in air --- \n");
	
	peakEc = 0.0;
	peakPc = 0.0;	
	for (int i=0;i<nData;i++) {
	        // printf ("%e %e %e %e %e %e\n",z[i],a[1],a[2],a[3],Emin,a[4]); 
		Dc[i] = P(z[i],a[1],a[2],a[3],Emin,a[4]);
		// printf ("Dc[i] = %e\n", Dc[i]);
		if(Dc[i] >= peakPc) {
			 peakPc = Dc[i];
			 peakEc = z[i];
		}
	}
			
	for (int i=0;i<nData;i++) {
		Dc[i] = Dc[i]/peakPc;
		if (dFlag != 1) printf ("%e %e %e\n", z[i], Dm[i], Dc[i]); 		
	}
	
	
   // free memory
   free_matrix(alpha,1,MA,1,MA);
   free_matrix(covar,1,MA,1,MA);
   free_vector(sig,1,NPT);
   free_vector(y,1,NPT);
   free_vector(x,1,NPT);
   free_vector(a,1,MA);
   free_ivector(ia,1,MA);
				
	return(PASS);
}

