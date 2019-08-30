/*************************************************************************
  Program:     model.c
  --------
  Description:
  ----------- 
	Simulate carbon and water cycles for a single pixel on a day.
***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:   J. Liu      
  Modified by:  X.F. Feng
  Last update:  July 2003
*****************************************************************************/
   

#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 

    
__device__  void gpumodel(float CNl,long jday,int long pix,long index,float lat_p,float lai_p,int short  lc_p,float CI,float TI,float sw_o1,float sw_o2,float sw_o3,
                           struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[],  float  x[],  float  b[], short climatedata[], 
						   struct     xvalue xx[], float co,float*soilw1_new,float* soilw2_new,float* soilw3_new,float *tmean,float *t_d,float *t_n,float *d_l)
    
	{
    
	//float *sw_n1,*sw_n2,*sw_n3;	
	float z[SIZEZ];
    float g[SIZEG];
    float  swn1,swn2,swn3;
    float  awc_p=0.5;
    float Tmean;
 //lixuansong
    if(jday==jday_start)
    {
    gpureadxx(pix,x,xx);	                	/* read snowpack and others */
	gpusetx(lc_p,awc_p,x,b);	                /* set x[2]-x[10] */ 		
    }

    /* read & compute forest-bgc daily climate data */ 
    else
    {
		gpureadxx(pix,x,xx);
    }

	gpuzcomp(jday,pix,lat_p,lc_p,CI,TI,awc_p,lai_p,b,x,z, climatedata,&Tmean); 

    *tmean=(z[4]+z[5])/2.0;

       *t_d=z[14];
	   *t_n=z[15];
	   *d_l=z[18];



	/* daily water and carbon dynamics */ 

	gpudoflux(CNl,pix,index,b,g,x,z,co,TI,sw_o1,sw_o2,sw_o3,HY1,HY2,HY3,&swn1,&swn2,&swn3); 


    (*soilw1_new)=swn1;(*soilw2_new)=swn2;(*soilw3_new)=swn3;


	gpuwritexx(pix,x,xx); 


/***********************************************************************/
//	free(g);
   	return;		
} 
 

  __device__  void gpusetx( int short lc_p, float awc_p, float x[], float b[])
    {  
    //x[2]=awc_p*1.1;  

/* Set x[3] to x[7], x[11]-x[16] to zero */ 

    x[3]=0.0; 
    x[4]=0.0; 
    x[5]=0.0; 
    x[6]=0.0; 
    x[7]=0.0; 
  
    x[11]=0.0; 
    x[12]=0.0;
	x[13]=0.0;
	x[14]=0.0;
	x[15]=0.0;
	x[16]=0.0;
	x[17]=0.0;
	x[18]=0.0;
    x[19]=0.0;

	 switch(lc_p)
    {
		 

//******* for MODIS, July ********* 
	 case 4: case 10:  //ENF
	x[9]=x[9]/b[1]*0.081;// 	x[9]=x[9]/b[1]*0.081;annavg livewood:annmax leaf(kgC/kgC) .0.081
        //if (x[9]==0) x[9]=30;	
   	x[10]=0.2317*x[9];					// root biomass in t/h b 

	break;
 case 1: case 7: case 8://EBF

		x[9]=x[9]/b[1]*0.162;// *0.162;annavg livewood:annmax leaf(kgC/kgC) 0.162  
		x[10]=exp(0.359)*pow(x[9],(float)0.639);

	break;
case 5: //DNF

		x[9]=x[9]/b[1]*0.152;// annavg livewood:annmax leaf(kgC/kgC) 0.152
		x[10]=0.2317*x[9];	
	
	break;
case 2: case 3://DBF
	
		x[9]=x[9]/b[1]*0.203;// annavg livewood:annmax leaf(kgC/kgC) 0.203
		x[10]=exp(0.359)*pow(x[9],(float)0.639);	
	
	break;

case 6: case 9: //MF
	
		x[9]=x[9]/b[1]*0.132;// annavg livewood:annmax leaf(kgC/kgC) 0.132
		x[10]=0.5*(0.2317*x[9]+exp(0.359)*pow(x[9],(float)0.639));	
	
	break;
 
case 11:case 12:case 14:  // 8:closed shrub;9:open shrub; 7: woody savanna;6:savannas;11: permanent wetlands
	
		x[9]=x[9]/b[1]*0.040;// annavg livewood:annmax leaf(kgC/kgC) 0.040
		x[10]=0.2317*x[9];
	
	break;

default: 
	
		x[9]=1.0;
		x[10]=0.2; // defult root biomass //zhoulei 4.1

    }
} 



