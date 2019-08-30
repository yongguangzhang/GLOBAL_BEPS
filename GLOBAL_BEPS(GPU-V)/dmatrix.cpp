/********************************************************************
**  Function Name: dmatrix
**  Description:   Allocates memory for float precision matrix
**
**  Details:       If the start position of the matrix is zero,
**                 a matrix of size (nrh,nch) is allocated, otherwise it
**                 allocates a dble matrix with range [nrl..nrh][ncl..nch]
**
**  References:
**                 Press, W.H., B.P. Flannery, S.A. Teukolsky, and
**                 W.T. Vetterling, "Numerical Recipes in C - The Art of
**                 of Scientific Computing," Cambridge University Press,
**                 1988.
*********************************************************************
**  Copyright:     CCRS
**  Written by:    David Fraser
**  Modified by:   X.F. Feng
**  Last update:   July 2003
******************************************************************/
#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 

float **dmatrix(

int *err,
int nrl, int nrh,int ncl,int nch)
{
        int i;
        float **m;

        *err = 0;	               /*NOERR*/

/*      Allocate pointers to rows */
        
		if (nrl ==0)    //当从第一行开始时，分配nrh个空间， nrh为总行数
          m = (float **) malloc((unsigned) (nrh)*sizeof(float));   

/*      Else allocate a dble matrix with range [nrl..nrh][ncl..nch] */
        else                 //当从nrh行开始时，分配(nrh-nrl+1)个空间
          m = (float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float));

        if (!m) {
	  printf("Error occured in allocating the matrix\n");   //// ？？当m=0时
          return (float **)0;
        }

        m -= nrl;    //  m=m-nrl（剩余行）  干嘛用

/*      Allocate rows and set pointers to them */
        if (ncl==0)   //如果是从第一列开始，则分配(nch)个空间，nch为总列数
          for (i=nrl; i<nrh; i++)
          {
            m[i] = (float *) malloc((unsigned) (nch)*sizeof(float));
            if (!m[i]) {
	      printf("Error occured in allocating the matrix\n");
              return (float **)0;
            }
            m[i] -= ncl;        //m[i]=m[i]-ncl    ??????????????什么用处（剩余列）
          }
        else
          for (i=nrl; i<=nrh; i++)
          {
            m[i] = (float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
            if (!m[i]) {
	      printf("Error occured in allocating the matrix\n");
              return (float **)0;
            }
            m[i] -= ncl;     //m[i]=m[i]-ncl 
          }

        return m;
}
