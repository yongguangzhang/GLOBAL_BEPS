/*************************************************************************
  Program:     readlaimax.c
  --------
  Description:
  -----------
  Read biomass from a file in a given line.

 ***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:  J. Liu       
  Modified by: X.F. Feng
  Last updatd:  July 2003
*****************************************************************************/

#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 

void readlaimax(int long line,struct xvalue xx[],int short *lc,int short *rcode,long tmpyear)
{
	FILE *in;                   /* Pointers biomass file */
	char filename[255];             ///改zfm3.26= "G:\\LAI_smooth\\00laimax.dat"
	//short int* laimaximg;   ///改   unsigned char
	 unsigned char* laimaximg; 
	long ptr;
	int size, nitems;
	int bytes=1;   ///改 short int		   /* Number of bytes for a unsigned char data type */
	int long pix;
    float omega;
//    float a1;


  laimaximg=(unsigned char *)malloc(npixels*sizeof(unsigned char)); ///改
 //laimaximg=(short int *)malloc(npixels*sizeof(short int)); ///改

  /*
 if(tmpyear<1980) sprintf(filename, "\\Global_Input\\Global_LAI_max\\%d%s",1980,"_S_LAI_max.dat"); //合成文件名
            else  sprintf(filename, "\\Global_Input\\Global_LAI_max\\%d%s",tmpyear,"_S_LAI_max.dat"); //合成文件名
*/
			
		if(tmpyear<1981) sprintf(filename, "\\Global_Input\\Global_LAI_max\\GlobMapLAIV3_%d%s",1982,"_LAImax.dat"); //合成文件名
                   else  sprintf(filename, "\\Global_Input\\Global_LAI_max\\GlobMapLAIV3_%d%s",tmpyear,"_LAImax.dat"); //合成文件名	
			
	// sprintf(filename, "\\Global_Input\\Global_LAI_max\\%d%s",1982,"_S_LAI_max.dat");//==================================================???

   if((in= fopen(filename,"rb")) ==NULL) {
     printf("Unable to open file %s\n", filename);
   *rcode = ERROR;
   return;
     }

   //Set pointer to start of line data 
     ptr=sizeof( char)*(line*npixels + pix_offset);  ///改
      
		fseek(in,ptr,0);     
        fread(&laimaximg[0],size=bytes,nitems=npixels,in);
		//fread(&laimaximg[0],sizeof(short),npixels,in);
/*
		if(line==9)
		{
			FILE* fout3=fopen("out33.txt","w");
			fprintf(fout3,"%d %d\n",ptr,sizeof(char));
			for(pix=pix_offset;pix<npixels;pix++){
				fprintf(fout3,"%d %d %d\n",pix,lc[pix],(int)laimaximg[pix]);
			}
			fclose(fout3);
		}
*/

  	for (pix=pix_offset; pix<npixels; pix++) { 


  		 /*
	switch(lc[pix])
    {
//************************** for MODIS, July***********************zhoulei 4.1
	case 6 : case 9://mixe =0.6, for MODIS 
		omega=0.65;           //0.6   //2013-0828由0.7改为0.65 
		break;

	case 1: case 2: case 3:case 7: case 8:    //broadleaf forest =0.7, for MODIS 
		omega=0.75;      //0.7       //2013-0828由0.8改为0.75    
		break;

	case 4: case 5: case 10: //coni=0.5, for MODIS 
		omega=0.55;   //0.55  7-21         //2013-0828由0.6改为0.55  
		break;


    case 11:  case 12: case 14:     //shrub=0.5, for MODIS, shrub same as conif, closed shrub and woody savanna
		omega=0.5;
		break;

	case 13: case 16: case 17: case 18: 	// grass=0.9, for MODIS, opened shrub & savanna & grassland
		omega=0.9;
		break;
	
	default:
		omega=0.8;
}  
	*/		
	xx[pix].x9=(float)laimaximg[pix]*0.1;//omega;  //3.25this is maximum value of LAI in one year
    xx[pix].x25=xx[pix].x9;
	} //end of pix loop
	
	

		fclose(in);
		free(laimaximg);

		return;
}






