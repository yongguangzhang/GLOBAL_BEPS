/*************************************************************************
  Program:     readclim.c
  --------
  Description:
  -----------
  read climate data. (grid file, interpolated  using ANUSPLIN, 687 weather stations in China. 


  Notes:
  Tmin: 16s, unit is 0.01C.
  Tmax: 16s, unit is 0.01C.
  Precip: 16s, Unit is 0.1mm.
  Rehum: 8u, unit is %.
  Radiation: 16s, unit is kj/m2.day

***************************************************************************
  LREIS (IGSNRR of CAS, China)
  Written by:   X. F. Feng       
  Last update:  July 2003
*****************************************************************************/

#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 
#include "malloc.h"


void readclim(  
long year, long jday,long line,	  
short climatedata[],
int short *rcode)
{
    
	char name_rh[255];                //文件名， Relative Humidity, daily, 8u bit
    char name_tmax[255];              //文件名， Maxmium Temperature, daily, 16s bit 
    char name_tmin[255];              //文件名，  Minmium Temperature, daily, 16s bit 
    //char name_pre[255];               //文件名，  Precipitaion, daily, 16s bit 
	char name_rad[255];               //文件名，  Total Solar Radiation, daily, 16u, bit 
   // char tmaxf[100],tminf[100],precipf[100],rehumf[100]; //sradf[100],
    FILE  *f_rh,*f_tmax,*f_tmin;// *f_pre,*f_rad ;  /* Climate data, Feng, June, 2003*/
   	long ptr;
   
  
//=========================================Prinston Climate data======================================

	    //合成文件名
	   // strcpy(tmaxf,   "\\Global_Input\\GlobalM_Prinston_2012_11\\M");
	  
		
		
		
       //strcpy(tmaxf,   "\\Global_Input\\GlobalM_Princeton_2015_07\\M");

strcpy(tmaxf,   "\\Global_Input\\GlobalM_CRU_2017_09V8\\M");

		if(jday<10){
	    sprintf(name_tmax,"%s%d%s%d%s%d%s",tmaxf, year,"\\Md", year,"00",jday,".dat");   ///zfm3.25
		}
	
		if(jday>=10 && jday<100) {
      sprintf(name_tmax,"%s%d%s%d%s%d%s",tmaxf, year,"\\Md", year,"0",jday,".dat");   ///zfm3.25
			}

		if(jday>=100){
	  sprintf(name_tmax,"%s%d%s%d%d%s",tmaxf, year,"\\Md", year,jday,".dat");   ///zfm3.25
	  	}
       //合成文件名结束

 

       if((f_tmax=fopen(name_tmax,"rb"))==NULL){
        printf("Can not open the file %s\n",name_tmax);    //从name_tmax文件中打开tmax
       exit(0);
	  }

ptr=((unsigned long)line*(unsigned long)npixels + (unsigned long)pix_offset);  //org   
fseek(f_tmax, 1L*ptr*5*sizeof(short),SEEK_SET);
fread(&climatedata[0],sizeof(short),npixels5,f_tmax); 

//if(year>1901) printf("year=%d jday=%d Tmax=%d Tmin=%d\n",year,jday,climatedata[43033],climatedata[47983]);
fclose(f_tmax);
 
//free(climatedata);
 

}
