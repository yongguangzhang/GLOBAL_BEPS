/*************************************************************************
  Program:     melt.c
  --------
  Description:
  ----------- 
	Compute snowmelt.
  Details:
  -------
	input:
	b[6]: snow melt temp. coefficient (m/C/day)(=0.001) 
	b[8]: snow aborptivity (=0.12)
	x[1]: snowpack (m)
	z[8]: shortwave radiation at forest floor (kj/m2/d)
	z[14]: daytime average temperature (C deg)		
	output:
	g[3]: total melt (m)
	g[8]: temperature melt (m)
	g[9]: radiation melt (m)
	g[13]: sublimation from snow (m)
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
}*/
__device__ void gpumelt (
    float b[],float g[],float x[],float z[])
    {
	float lh_fus=3.5E+5;
//	float lh_evp=2.5E+6; 
	float lh_sub=2.8E+6;   ///?????????3.25
	short int lc_p;
/* Snowmelt Temp Coeff is different for different land cover type  X.F. July */

	lc_p=(int short )z[23];
    switch(lc_p)
	{

	/* for LCTM, reference about Snowmelt Temp Coffe. is  G.W. Kite et al. 1988 X.F. July */

	case 4: case 10: case 5:              /* conifer X.F. Aug. */ // conifer ;zhoulei 4.1
		b[6]=0.001/RTIMES;		           /* unit m/degree C/d  */     //2018-03-05: 0.0022 to 0.0010
		break;

	case 1: case 7:case 8: case 2:case 3:            /* braodleaf forest,for LCTM */ // braodleaf forest zhoulei 4.1
        b[6]=0.0010/RTIMES;		           /* unit m/degree C/d  */ 
		break;

	case 6: case 9:                                /*mixe ,for LCTM */ //mixed zhoulei 4.1
		b[6]=0.0015/RTIMES;		           /* unit m/degree C/d  */ 
		break;
		
	case 11:case 12:case 14:  /*closedshrub, pasture and crop for LCTM */ //closed shrub &grassland &cropland ;zhoulei 4.1
		b[6]=0.0020/RTIMES;		           /* unit m/degree C/d  */ 
		break;

	case 13:  case 16:case 17: case 18: //zfm3.25  10,12,14
        b[6]=0.001/RTIMES;
		break;

	default:                             /* other land cover */
	b[6]=0.006/RTIMES;
	}



/*	b[8]: snow absorption */
	b[8]=0.3/3.0*0.1;


    if (z[14] > 0.0) { 
	    g[8] = z[8]*b[8]/lh_fus;	/* radiation  melt */ 
		g[9] = b[6] * z[14];		/* temperature melt */ 
		g[3] = __min (x[1], (g[8]+g[9]));
		g[13] =0.0;
	}
    else { 
		g[8]=0.0;
		g[9]=0.0;
		g[3]=0.0;
 	
		g[13] = (z[8]-z[19])*b[8]/lh_sub*0.001;  //2012-12-04
       	//	g[13]=z[19]*b[8]/lh_sub;
        g[13]=__max(0,g[13]);        
		g[13] = __min(x[1],g[13] );
		

	}


	return;

	}


