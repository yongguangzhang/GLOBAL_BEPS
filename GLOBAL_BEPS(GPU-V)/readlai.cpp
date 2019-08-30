/*************************************************************************
  Function:     readlai.c
  --------
  Description:
  -----------
 	Read LAI files
***************************************************************************
  CCRS (EMS/Applications Division) 
  Written by:   J. Liu
  Modified by: X.F. Feng
  Last updatd:  July 2003
*****************************************************************************/

#include"beps.h"

void readlai(int long line, float *lai,long jday, int short *rcode,long tmpyear)
{
	FILE *infile,*infile1,*infile2,*infile3,*infile4,*infile5,*infile6,*infile7,*infile8,*infile9,*infile10;
	int long pix;		
  	int long ptr;
  	int long size;
  	int long nitems;
  	char filename[255];
  
    char filename1[255];
    char filename2[255];
    char filename3[255];

    char filename4[255];
    char filename5[255];

    char filename6[255];
    char filename7[255];
    char filename8[255];

    char filename9[255];
    char filename10[255];
     float *laiimg0;
    unsigned char * laiimg;   /// 改 unsigned char
	 
	unsigned char * laiimg1;   /// 改 unsigned char
	 
	unsigned char * laiimg2;   /// 改 unsigned char
	   
	unsigned char * laiimg3;   /// 改 unsigned char
	 unsigned char * laiimg4;   /// 改 unsigned char
	   
	unsigned char * laiimg5;   /// 改 unsigned char
	
    unsigned char * laiimg6;   /// 改 unsigned char
	 
	unsigned char * laiimg7;   /// 改 unsigned char
	 
	unsigned char * laiimg8;   /// 改 unsigned char
	   
	unsigned char * laiimg9;   /// 改 unsigned char
	unsigned char * laiimg10;   /// 改 unsigned char
	   
	int JDY;
    int YEAR;
 



   JDY=jday;

   if(tmpyear<=2000) {	


	   if(tmpyear<1981){   //2013-11-11 由1981改为1982
		  
		laiimg1=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改
		laiimg2=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改
		laiimg3=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改
       laiimg4=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改
	laiimg5=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改

        laiimg6=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改
		laiimg7=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改
		laiimg8=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改
        laiimg9=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改
		laiimg10=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改

YEAR=1981;
if(YEAR%4==0 && jday>=60) JDY=jday+1;
else  JDY=jday;
		   if(JDY<10)             sprintf(filename1, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
		   if(JDY>=10 && JDY<100) sprintf(filename1, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
		   if(JDY>=100)           sprintf(filename1, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");

		   YEAR=1982;

		   if(YEAR%4==0 && jday>=60) JDY=jday+1;
else  JDY=jday;

if(JDY<10)             sprintf(filename2, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
if(JDY>=10 && JDY<100) sprintf(filename2, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
if(JDY>=100)           sprintf(filename2, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");


		   YEAR=1983;
		   if(YEAR%4==0 && jday>=60) JDY=jday+1;
else  JDY=jday;

if(JDY<10)             sprintf(filename3, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
if(JDY>=10 && JDY<100) sprintf(filename3, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
if(JDY>=100)           sprintf(filename3, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");
          
		   YEAR=1984;
		   if(YEAR%4==0 && jday>=60) JDY=jday+1;
else  JDY=jday;

/*
if(JDY<10)             sprintf(filename4, "\\Global_Input\\Global_LAI\\2017LAI_smoothed\\GlobMapLAIV3_S_C20180128_%d%s%d%s",YEAR,"00",JDY,".dat");
if(JDY>=10 && JDY<100) sprintf(filename4, "\\Global_Input\\Global_LAI\\2017LAI_smoothed\\GlobMapLAIV3_S_C20180128_%d%s%d%s",YEAR,"0",JDY,".dat");   
if(JDY>=100)           sprintf(filename4, "\\Global_Input\\Global_LAI\\2017LAI_smoothed\\GlobMapLAIV3_S_C20180128_%d%d%s",YEAR,JDY,".dat");
*/
if(JDY<10)             sprintf(filename4, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
if(JDY>=10 && JDY<100) sprintf(filename4, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
if(JDY>=100)           sprintf(filename4, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");


		   YEAR=1985;
		   if(YEAR%4==0 && jday>=60) JDY=jday+1;
		 else  JDY=jday;
		 if(JDY<10)             sprintf(filename5, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
		 if(JDY>=10 && JDY<100) sprintf(filename5, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
		 if(JDY>=100)           sprintf(filename5, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");

		   YEAR=1986;
		   if(YEAR%4==0 && jday>=60) JDY=jday+1;
		 else  JDY=jday;
		 if(JDY<10)             sprintf(filename6, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
		 if(JDY>=10 && JDY<100) sprintf(filename6, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
		 if(JDY>=100)           sprintf(filename6, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");;



		   YEAR=1987;
		   if(YEAR%4==0 && jday>=60) JDY=jday+1;
		  else  JDY=jday;
		  if(JDY<10)             sprintf(filename7, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
		  if(JDY>=10 && JDY<100) sprintf(filename7, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
		  if(JDY>=100)           sprintf(filename7, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");

		   YEAR=1988;
		   if(YEAR%4==0 && jday>=60) JDY=jday+1;
		   else  JDY=jday;
		   if(JDY<10)             sprintf(filename8, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
		   if(JDY>=10 && JDY<100) sprintf(filename8, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
		   if(JDY>=100)           sprintf(filename8, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");



		   YEAR=1989;
		   if(YEAR%4==0 && jday>=60) JDY=jday+1;
		  else  JDY=jday;
		  if(JDY<10)             sprintf(filename9, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
		  if(JDY>=10 && JDY<100) sprintf(filename9, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
		  if(JDY>=100)           sprintf(filename9, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");
		  
		  
		  
		  YEAR=1990;
		   if(YEAR%4==0 && jday>=60) JDY=jday+1;
		   if(JDY<10)             sprintf(filename10, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"00",JDY,".dat");
		   if(JDY>=10 && JDY<100) sprintf(filename10, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",YEAR,"0",JDY,".dat");   
		   if(JDY>=100)           sprintf(filename10, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",YEAR,JDY,".dat");
   
	   //Open LAI file 
	   if ((infile1=fopen(filename1, "rb")) == NULL) {
		   printf("Unable to open file %s\n", filename1);
		   *rcode = ERROR;
		   return;
	   }


	   if ((infile2=fopen(filename2, "rb")) == NULL) {
		    printf("Unable to open file %s\n", filename2);
		   *rcode = ERROR;
		   return;
	   }

 
	   if ((infile3=fopen(filename3, "rb")) == NULL) {
		     printf("Unable to open file %s\n", filename3);
		   *rcode = ERROR;
		   return;
	   }


	   if ((infile4=fopen(filename4, "rb")) == NULL) {
		   printf("Unable to open file %s\n", filename4);
		   *rcode = ERROR;
		   return;
	   }
 
	   if ((infile5=fopen(filename5, "rb")) == NULL) {
		     printf("Unable to open file %s\n", filename5);
		   *rcode = ERROR;
		   return;
	   }


        if ((infile6=fopen(filename6, "rb")) == NULL) {
		    printf("Unable to open file %s\n", filename6);
		   *rcode = ERROR;
		   return;
	   }


	   if ((infile7=fopen(filename7, "rb")) == NULL) {
		     printf("Unable to open file %s\n", filename7);
		   *rcode = ERROR;
		   return;
	   }

 
	   if ((infile8=fopen(filename8, "rb")) == NULL) {
		  printf("Unable to open file %s\n", filename8);
		   *rcode = ERROR;
		   return;
	   }

	 
	   if ((infile9=fopen(filename9, "rb")) == NULL) {
		    printf("Unable to open file %s\n", filename9);;
		   *rcode = ERROR;
		   return;
	   }
 
	   if ((infile10=fopen(filename10, "rb")) == NULL) {
		  printf("Unable to open file %s\n", filename10);
		   *rcode = ERROR;
		   return;
	   }
  
	ptr=(line*npixels + pix_offset)*sizeof(char);   ///改
	fseek(infile1,ptr,0);
	fread(&laiimg1[0],size=1,nitems=npixels,infile1);  ///改size=1
	
	fseek(infile2,ptr,0);
	fread(&laiimg2[0],size=1,nitems=npixels,infile2);  ///改size=1
	
	fseek(infile3,ptr,0);
	fread(&laiimg3[0],size=1,nitems=npixels,infile3);  ///改size=1
	

    fseek(infile4,ptr,0);
	fread(&laiimg4[0],size=1,nitems=npixels,infile4);  ///改size=1
	
    fseek(infile5,ptr,0);
	fread(&laiimg5[0],size=1,nitems=npixels,infile5);  ///改size=1
	

    fseek(infile6,ptr,0);
	fread(&laiimg6[0],size=1,nitems=npixels,infile6);  ///改size=1
	
	fseek(infile7,ptr,0);
	fread(&laiimg7[0],size=1,nitems=npixels,infile7);  ///改size=1
	
	fseek(infile8,ptr,0);
	fread(&laiimg8[0],size=1,nitems=npixels,infile8);  ///改size=1
	

    fseek(infile9,ptr,0);
	fread(&laiimg9[0],size=1,nitems=npixels,infile9);  ///改size=1
	
    fseek(infile10,ptr,0);
	fread(&laiimg10[0],size=1,nitems=npixels,infile10);  ///改size=1


	for (pix=pix_offset;pix<npixels;pix++)
	{
		
		//if(laiimg2[pix]<150 && laiimg2[pix]>0) {
		lai[pix]=laiimg1[pix]+laiimg2[pix]+laiimg3[pix]+laiimg4[pix]+laiimg5[pix]+laiimg6[pix]+laiimg7[pix]+laiimg8[pix]+laiimg9[pix]+laiimg10[pix];
		lai[pix]=laiimg2[pix]+laiimg3[pix]+laiimg4[pix]+laiimg5[pix]+laiimg6[pix];//+laiimg7[pix]+laiimg8[pix]+laiimg9[pix]+laiimg10[pix];
		
		lai[pix]=lai[pix]/50.0;


		//}
		//else                 lai[pix]=0.05;                 //居20111-4-23改;
		
	
		if(lai[pix]<0.01) lai[pix]=0.01;
		if(lai[pix]>10.0)  lai[pix]=10.0;

	}
	   
   	   fclose(infile1); 	   fclose(infile2);   	   fclose(infile3);
	  
	   
	  fclose(infile4);
	   fclose(infile5);

       fclose(infile6);
	   fclose(infile7);
	   fclose(infile8);
	   fclose(infile9);
	   fclose(infile10);

	   free(laiimg1); free(laiimg2); free(laiimg3);

	  
	free(laiimg4); free(laiimg5); free(laiimg6); free(laiimg7); free(laiimg8);free(laiimg9); free(laiimg10);
	
 }
	   else{   //=========================================================================Years from 1981 to 1999==============================================
       
		              laiimg=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改
 
					  if(tmpyear%4==0 && jday>=60) JDY=jday+1;

	   if(JDY<10)      sprintf(filename, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",tmpyear,"00",JDY,".dat");  
if(JDY>=10 && JDY<100) sprintf(filename, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%s%d%s",tmpyear,"0",JDY,".dat"); 
if(JDY>=100)           sprintf(filename, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_C20180222_%d%d%s",tmpyear,JDY,".dat"); 
 

	  
	   if ((infile=fopen(filename, "rb")) == NULL) { //===========================================================Open LAI file===============================
		  printf("Unable to open file %s\n", filename);
		   *rcode = ERROR;
		   return;
	   }

	   
	ptr=(line*npixels + pix_offset)*sizeof(char);   ///改
	fseek(infile,ptr,0);
	fread(&laiimg[0],size=1,nitems=npixels,infile);  ///改size=1
	
	
	for (pix=pix_offset;pix<npixels;pix++)
	{
		
		if(laiimg[pix]<150) {
	    lai[pix]=laiimg[pix]*0.1;  
	   //  if(tmpyear==1991||tmpyear==1992)lai[pix]=lai[pix]*1.02;	
		}
		
		else    lai[pix]=0.01;                 //居20111-4-23改;
		
		if (laiimg[pix]==0) lai[pix]=0.01;     //2009-1-1-22
		
		if(lai[pix]<0.01) lai[pix]=0.01;
		if(lai[pix]>10.0)  lai[pix]=10.0;
	}
	   

 	fclose(infile);
   free(laiimg);
      }
  
   }
else{   //===========================================Years after 2000

laiimg=(unsigned char*)malloc(npixels*sizeof(unsigned char)); ///改

	if(JDY<10)             sprintf(filename, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_%d%s%d%s",tmpyear,"00",JDY,".dat");  
	if(JDY>=10 && JDY<100) sprintf(filename, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_%d%s%d%s",tmpyear,"0",JDY,".dat"); 
	if(JDY>=100)           sprintf(filename, "\\Global_Input\\Global_LAI\\20180222\\GlobMapLAIV3_S_%d%d%s",tmpyear,JDY,".dat"); 


	if ((infile=fopen(filename, "rb")) == NULL) {
	  printf("Unable to open file %s\n", filename);
		*rcode = ERROR;
		return;
	}

	ptr=(line*npixels + pix_offset)*sizeof(char);   ///改
	fseek(infile,ptr,0);
	fread(&laiimg[0],size=1,nitems=npixels,infile);  ///改size=1


	for (pix=pix_offset;pix<npixels;pix++)
	{
		if(laiimg[pix]<150)    lai[pix]=laiimg[pix]*0.1;          ///zfm09-2-2, (laiimg[pix]-1)*0.01, 3.25
		else                   lai[pix]=0.01;                 //居20111-4-23改;
		if (laiimg[pix]==0)    lai[pix]=0.01;                 //2009-1-1-22
		 
		if(lai[pix]<0.01) lai[pix]=0.01;
		if(lai[pix]>10.0)  lai[pix]=10.0;
       
	   //if(tmpyear==2012)  lai[pix]=lai[pix]*1.04;   //1.08
      // if(tmpyear==2010) lai[pix]=lai[pix]*0.99;  //0.95
      // if(tmpyear==2006 || tmpyear==2009 )lai[pix]=lai[pix]*1.01;  //1.05
		}
	//printf("tmpyear=%d %d %f\n",tmpyear,JDY,lai[2540]);
	fclose(infile);
	free(laiimg);
}



}
