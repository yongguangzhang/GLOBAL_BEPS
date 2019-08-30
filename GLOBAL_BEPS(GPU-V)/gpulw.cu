/*************************************************************************
  Program:     doflux.c
  --------
  Description:
  ----------- 
	Calculate and update the fluxes in carbon and water cycles.
***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:   J. Liu      
  Last update:	December 2000
*****************************************************************************/
	 
#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 


   __device__ void gpulw1(float FC1,float FC2, float FC3, float sw1, float  sw2, float sw3, float *lw1,float *lw2, float *lw3)
    {          
 
float m,Msat, Mopt,B1,B2,B3,B; 


B3=100.0;
//Layer1
m=-0.0008*FC1*FC1+0.0603*FC1-0.6769+0.5;

Mopt=-0.0062*FC1*FC1+1.0745*FC1+26.515;

Msat=(-0.00006*FC1*FC1+0.0124*FC1+0.0994);

B1=pow(sw1,m)-pow(Mopt,m);
B2=pow(Mopt,m)-pow(B3, m);
B=(B1/B2)*(B1/B2);

*lw1=0.75*pow(Msat,B)+0.25;

//Layer2
m=-0.0008*FC2*FC2+0.0603*FC2-0.6769+0.5;

Mopt=-0.0062*FC2*FC2+1.0745*FC2+26.515;
Msat=(-0.00006*FC2*FC2+0.0124*FC2+0.0994);

B1=pow(sw2,m)-pow(Mopt,m);
B2=pow(Mopt,m)-pow(B3, m);
B=(B1/B2)*(B1/B2);

*lw2=0.75*pow(Msat,B)+0.25;


//Layer3
m=-0.0008*FC3*FC3+0.0603*FC3-0.6769+0.5;

Mopt=-0.0062*FC3*FC3+1.0745*FC3+26.515;
Msat=(-0.00006*FC3*FC3+0.0124*FC3+0.0994);

B1=pow(sw3,m)-pow(Mopt,m);
B2=pow(Mopt,m)-pow(B3, m);
B=(B1/B2)*(B1/B2);

*lw3=0.75*pow(Msat,B)+0.25;

    return;
}





