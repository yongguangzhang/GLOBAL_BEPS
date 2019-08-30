/*************************************************************************
  Program:     inter.c
  --------
  Description:
  ----------- 
	Compute rain, snow to ground and to air.
  Details:
  -------
	input:
	b[4]: snow & rain interception coefficient (=0.00025 m/lai/day)
	z[3]: precipitation
	z[8]: daily total short wave radiation (kj/m2/d)，=J/(3.6*m2*s)=(w/m2)/3.6??????????
	z[10]: lai
	z[14]: daytime average temperature (C deg)		
	output:
	g[1]:ground rain 
	g[2]:ground snow
	g[4]:evaporation of precipitation from canopy (m of h2o/ha/day)	
	g[5]:sublimation of precipitation from canopy (m of h2o/ha/day)		
	g[6]:ground rain from canopy
	g[7]:ground snow from canopy
	g[11]:intercepted precip. 
***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:   J. Liu      
  Modified by: X.F. Feng
  Last update:  July 2003
*****************************************************************************/
	
#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 

/*float min(float a,float b)
{
	if (a<b)
		return a;
	return b;
}
float max(float a,float b)
{
	if (a>b)
		return a;
	return b;
}
  */
__device__ void gpuinter(
    float  b[],float g[],float x[],float z[])
{
    // float lh_evp=2.5E+6; 
    // float lh_sub=2.8E+6;
    float  Eva_inter_sunlit,Eva_inter_shaded, Eva_inter,Rn;
float  a1;

	/* fresh snow */
	
	b[8]=0.05; 
 
/* intercepted precip */
	a1=b[4]*(1-exp(-0.5*z[10]))*5;     //*0.5;    // 2012-10-15日恢复, 并加系数0.5   z[10]为叶面积指数 
    g[11]= min(z[3], a1);              // 2012-10-15日恢复  z[3]为降水量    g[11] 冠层截留降水量 

//g[11]= __min(z[3],z[10]* b[4]*0.25);      // 2012-10-15日前用的方程



   
/* rain vs snow */     //g[1]为 
 
    if (z[14]  > 0.0 ) {    
        g[1] = z[3] - g[11];
		g[2] = 0.0; 
    } 
    else {
	    g[1]=0.0; 
        g[2] = z[3] - g[11];		
     }

/* water loss from canopy */
/*	
if (z[14]  > 0.0 ) {
        g[4] = min(g[11], z[8]*b[7]/lh_evp);
		g[5] =0.0;
 
    } 
    else {
		g[4] =0.0;
        g[5] = min(g[11], z[8]*b[8]/lh_sub);	
     }
*/

   


 
if (z[14]  > 0.0 ) {
    Rn=z[38]*z[32]+z[39]*z[33];

	Eva_inter_sunlit= z[18]*gpupenmon(z[14],z[16],Rn,g[40],0.05)*z[32];

    //Eva_inter_sunlit= z[18]*penmon(z[14],z[16],z[38],g[40],0.005)*z[32];

	Eva_inter_shaded= 0;//z[18]*penmon(z[14],z[16],z[39],g[40],0.005)*z[33];  

    Eva_inter=Eva_inter_sunlit+Eva_inter_shaded;
   
	if(g[11]>0) {
		g[90]=__min(1.0,g[11]/Eva_inter);     
	}
		else g[90]=0.0;
	
	   g[4] = __min(g[11], Eva_inter);   //2009_8_03 Ju 改动
	// g[4] = min(g[11], z[8]*b[7]/lh_evp); //原来的算法
	
	   
	   g[5] =0.0;
 
    } 
    else {
		g[4] =0.0;
    Rn=z[38]*z[32]+z[39]*z[33];

	Eva_inter_sunlit= z[18]*gpupenmon(z[14],z[16],Rn,g[40],0.05)*z[32];

    //Eva_inter_sunlit= z[18]*penmon(z[14],z[16],z[38],g[40],0.005)*z[32];

	Eva_inter_shaded= 0;//z[18]*penmon(z[14],z[16],z[39],g[40],0.005)*z[33];   

    Eva_inter=Eva_inter_sunlit+Eva_inter_shaded;

    if(g[11]>0) g[90]=__min(1.0,g[11]/Eva_inter);

	else g[90]=0.0;

	g[5] = __min(g[11], Eva_inter);   //2009_8_03 Ju 改动 

	//g[5] = __max(g[11], 0)	;


  if(g[5]<0.0)  g[5]=0;



		
	//	g[5] = min(g[11], z[8]*b[8]/lh_sub);//原来的算法	
     }



/* canopy water to ground */
	if (z[14]  > 0.0 ) {      
		g[6] = __max(0.0, g[11]-g[4]);
		g[7] =0.0;
    } 
    else { 
		g[6]=0.0;
		g[7] = __max(0.0, g[11]-g[5]);
     }

    return;
}
 
