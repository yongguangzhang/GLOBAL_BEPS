/*************************************************************************
  Program:     readconf.c
  --------
  Description:
  -----------
  Reads control file for NPP calculations.  The control file
  defines the area in which NPP will be calulated and the data
  files used in the calculations.

***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:	D. Fraser/J. Liu  
  Modified by: X.F. Feng
  Last update:  July 2003
*****************************************************************************/

#include "beps.h"

void readconf(char *conf,int *rcode)

{

	char field_descr[255];	   /* Field description in control file */
	FILE *fcon;		           /* Control file pointer */

	 
 
/*      If can not open control parameter file */
        if ((fcon = fopen(conf,"r")) == NULL) {
          printf("Can not open parameter file <%s>, exitting program ...\n",
                                                                conf);
          *rcode = ERROR;
	  return;
        }

/*	Read control file parameters */
    field_descr[0]='\0';                     // ????????????把字符串结束标志符为何放这里
   // fscanf(fcon,"%s %d", field_descr, &jday_start);    //从指针文件fcon读取jday_start
  //  field_descr[0]='\0';
 //   fscanf(fcon,"%s %d", field_descr, &jday_end);  zfm3.26
    field_descr[0]='\0';
	fscanf(fcon,"%s %d", field_descr, &cycle_display);
    field_descr[0]='\0';
	fscanf(fcon,"%s %f", field_descr, &factor);
    field_descr[0]='\0';
	fscanf(fcon,"%s %d", field_descr, &lin_offset);
    field_descr[0]='\0';
//	fscanf(fcon,"%s %d", field_descr, &nlines);
    field_descr[0]='\0';
	fscanf(fcon,"%s %d", field_descr, &pix_offset);
    field_descr[0]='\0';
//	fscanf(fcon,"%s %d", field_descr, &npixels);
	field_descr[0]=awcf[0]='\0';
	fscanf(fcon,"%s %s", field_descr, awcf);
	field_descr[0]=landcoverf[0]='\0';
	fscanf(fcon,"%s %s", field_descr, landcoverf);
	field_descr[0]=leafareaf[0]='\0';
	fscanf(fcon,"%s %s", field_descr, leafareaf);
    field_descr[0]=latitudef[0]='\0';
	fscanf(fcon,"%s %s", field_descr, latitudef);
	field_descr[0]=longitudef[0]='\0';
	fscanf(fcon,"%s %s", field_descr, longitudef);
/*    field_descr[0]=tmaxf[0]='\0';
	fscanf(fcon,"%s %s", field_descr, tmaxf);
    field_descr[0]=tminf[0]='\0';
	fscanf(fcon,"%s %s", field_descr, tminf);
	field_descr[0]=precipf[0]='\0';
	fscanf(fcon,"%s %s", field_descr, precipf);
	field_descr[0]=rehumf[0]='\0';
	fscanf(fcon,"%s %s", field_descr, rehumf);
	field_descr[0]=sradf[0]='\0';
	fscanf(fcon,"%s %s", field_descr, sradf);*/
	field_descr[0]=snowpackf[0]='\0';   
	fscanf(fcon,"%s %s", field_descr, snowpackf);
	field_descr[0]=soilwaterf[0]='\0';
	fscanf(fcon,"%s %s", field_descr, soilwaterf);
	field_descr[0]='\0'; 
	fscanf(fcon,"%s %c", field_descr, &output_NPP);	
	field_descr[0]='\0'; 
	fscanf(fcon,"%s %c", field_descr, &output_res);
	field_descr[0]='\0';
	fscanf(fcon,"%s %c", field_descr, &output_eva);	
	field_descr[0]='\0';
	fscanf(fcon,"%s %c", field_descr, &output_tra);
	field_descr[0]='\0';
	fscanf(fcon,"%s %c", field_descr, &output_site);

	
/*	Output control parameters */ //zhoulei 4.1
printf("\n Control Parameters\n\n");
printf(" Julian day_start [1,366]: %d\n",jday_start);
//printf(" Julian day_end [1,366]: %d\n",jday_end);
printf(" Cycle display [1,366]: %d\n",cycle_display);
printf(" NPP multiplication factor [20,20000]: %f\n",factor);
printf(" Line offset (PCI coordintate system) [0,429]: %d\n",lin_offset);
printf(" Number of Lines (PCI coordintate system) [1,430]: %d\n",nlines);
printf(" Pixel offset (PCI coordintate system) [0,529]: %d\n",pix_offset);
printf(" Number of pixels (PCI coordintate system) [1,530]: %d\n",npixels);
printf(" Available soil water holding capacity data file:  %s\n",awcf);
printf(" Land cover data file:  %s\n",landcoverf);
printf(" LAI data file:  %s\n",leafareaf);
printf(" Latitude data file:  %s\n",latitudef);
printf(" Longitude data file:  %s\n",longitudef);
/*printf(" Air maxmum temperature  data file:  %s\n",tmaxf);
printf(" Air minmum temperature daily data file:  %s\n",tminf);
printf(" Precipitation data file:  %s\n",precipf);
printf(" Relative humidity data file:  %s\n",rehumf);
printf(" Solair radiation data file:  %s\n",sradf);*/
printf(" Snowpack data file:  %s\n",snowpackf);
printf(" Soilwater data file:  %s\n",soilwaterf);  
printf(" Output NPP?  %c\n",output_NPP);
printf(" Output autorespiration?  %c\n",output_res);
printf(" Output evapiration?  %c\n",output_eva);
printf(" Output transpiration?  %c\n",output_tra); 
printf(" Output site file?  %c\n",output_site); 



fclose(fcon);


}



