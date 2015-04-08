#include <math.h>
#define NRANSI
#include "nrutil.h"
#ifdef NORM
#include <stdio.h>
#include <stdlib.h>
#endif

void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int), float *alamda)
{
	void covsrt(float **covar, int ma, int ia[], int mfit);
	void gaussj(float **a, int n, float **b, int m);
	void mrqcof(float x[], float y[], float sig[], int ndata, float a[],
		int ia[], int ma, float **alpha, float beta[], float *chisq,
		void (*funcs)(float, float [], float *, float [], int));
	int j,k,l,m;
	static int mfit;
	static float ochisq,*atry,*beta,*da,**oneda;

	if (*alamda < 0.0) {
		atry=vector(1,ma);
		beta=vector(1,ma);
		da=vector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=matrix(1,mfit,1,1);
		*alamda=0.001;
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			for (j++,k=0,m=1;m<=ma;m++) {
				if (ia[m]) {
					k++;
					covar[j][k]=alpha[j][k];
				}
			}
			covar[j][j]=alpha[j][j]*(1.0+(*alamda));
			oneda[j][1]=beta[j];
		}
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,ia,mfit);
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				for (j++,k=0,m=1;m<=ma;m++) {
					if (ia[m]) {
						k++;
						alpha[j][k]=covar[j][k];
					}
				}
				beta[j]=da[j];
				a[l]=atry[l];
			}
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
void covsrt(float **covar, int ma, int ia[], int mfit)
{
        int i,j,k;
        float swap;

        for (i=mfit+1;i<=ma;i++)
                for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
        k=mfit;
        for (j=ma;j>=1;j--) {
                if (ia[j]) {
                        for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
                        for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
                        k--;
                }
        }
}
#undef SWAP

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void gaussj(float **a, int n, float **b, int m)
{
        int *indxc,*indxr,*ipiv;
        int i,icol,irow,j,k,l,ll;
        float big,dum,pivinv,temp;

        indxc=ivector(1,n);
        indxr=ivector(1,n);
        ipiv=ivector(1,n);
        for (j=1;j<=n;j++) ipiv[j]=0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if (ipiv[j] != 1)
                                for (k=1;k<=n;k++) {
                                        if (ipiv[k] == 0) {
                                                if (fabs(a[j][k]) >= big) {
                                                        big=fabs(a[j][k]);
                                                        irow=j;
                                                        icol=k;
                                                }
                                        } else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
                                }
                ++(ipiv[icol]);
                if (irow != icol) {
                        for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
                        for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
                }
                indxr[i]=irow;
                indxc[i]=icol;
                if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
                pivinv=1.0/a[icol][icol];
                a[icol][icol]=1.0;
                for (l=1;l<=n;l++) a[icol][l] *= pivinv;
                for (l=1;l<=m;l++) b[icol][l] *= pivinv;
                for (ll=1;ll<=n;ll++)
                        if (ll != icol) {
                                dum=a[ll][icol];
                                a[ll][icol]=0.0;
                                for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
                                for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
                        }                                                               }
        for (l=n;l>=1;l--) {
                if (indxr[l] != indxc[l])
                        for (k=1;k<=n;k++)
                                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
        }
        free_ivector(ipiv,1,n);
        free_ivector(indxr,1,n);
        free_ivector(indxc,1,n);
}
#undef SWAP

void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[],
        int ma, float **alpha, float beta[], float *chisq,
        void (*funcs)(float, float [], float *, float [], int))
{
        int i,j,k,l,m,mfit=0;
        float ymod,wt,sig2i,dy,*dyda;

        dyda=vector(1,ma);
        for (j=1;j<=ma;j++)
                if (ia[j]) mfit++;
        for (j=1;j<=mfit;j++) {
                for (k=1;k<=j;k++) alpha[j][k]=0.0;
                beta[j]=0.0;
        }
        *chisq=0.0;
#ifdef NORM
		  float *yMod =(float *)malloc((ndata+1)*sizeof(float));
		  float yModMax = 0.0;
		  for (i=1;i<=ndata;i++) {
		  		(*funcs)(x[i],a,&ymod,dyda,ma);
				yMod[i] = ymod;
				// printf ("%f\t%f\n", x[i], yMod[i]); 
				if (ymod >= yModMax) yModMax = ymod;
		  }
		  /* Normalization to 100% */
		  if (yModMax <= 0.0) {
		  		printf ("ERROR: Normalization Factor %f <= 0.0\n", yModMax);
				exit(1);
		  }
		  float normScale = 100.0/yModMax; 
   	  for (i=1;i<=ndata;i++) {
		  		yMod[i] *= normScale;		 
		  }
#endif 		  
        for (i=1;i<=ndata;i++) {
                (*funcs)(x[i],a,&ymod,dyda,ma);
#ifdef NORM
					 ymod = yMod[i];
#endif
                sig2i=1.0/(sig[i]*sig[i]);
                dy=y[i]-ymod;
                for (j=0,l=1;l<=ma;l++) {
                        if (ia[l]) {
                                wt=dyda[l]*sig2i;
                                for (j++,k=0,m=1;m<=l;m++)
                                        if (ia[m]) alpha[j][++k] += wt*dyda[m];
                                beta[j] += dy*wt;
                        }
                }
                *chisq += dy*dy*sig2i;
        }
        for (j=2;j<=mfit;j++)
                for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
        free_vector(dyda,1,ma);
		  
#ifdef NORM
		  free(yMod);
#endif
}

#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5+4. */
