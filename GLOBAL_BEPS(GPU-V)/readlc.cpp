/*************************************************************************
  Function:     readlc.c 
  --------
  Description:
  -----------
  Reads landcover data for a specific line,  pixel dimensions
  Land cover data is stored in binary format.
***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:   J. Liu/D. Fraser
  Modified by: X.F. Feng
  Last updatd:  July 2003
*****************************************************************************/

#include"beps.h"

void readlc( int long line, int  short*lc,int short *rcode)
{		
	FILE *infile;
   int long  pix;		
  	int long  ptr;
  	//long int size, nitems;  
    unsigned char *lcimg;   /// 改 unsigned char
	//extern char landcoverf[255];
    // long aa;
 
	lcimg=(unsigned char *)malloc(npixels*sizeof(unsigned char)); 
/*	Open Land cover data file */

	//for (pix=0;pix<npixels;pix++) lcimg[pix]=0;
		
	if ((infile=fopen("\\Global_Input\\Global_cover\\Global_0727_landcover.raw", "rb")) == NULL) {	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", landcoverf);
	  *rcode = ERROR;
	  return;
    	}

/*	Calculate starting point of data */
	ptr=(line*npixels )*sizeof(char);   

   /*将文件指针移动到文件的开始*/
    fseek(infile,ptr,SEEK_SET);
	fread(&lcimg[0],1,npixels,infile);  

for (pix=0;pix<npixels;pix++) {
/*	Read data */
  	lc[pix]=lcimg[pix];
		}
 
/*	Close Land cover file */
free(lcimg);
  	fclose(infile);

	pix=1;
return;
}
