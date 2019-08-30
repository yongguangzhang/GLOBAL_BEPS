/*************************************************************************
  Program:     readindex.c
  --------
  Description:
  -----------
  Interpolates ncar (gridded) data using bilinear interpolation
  for all pixels in a line.

  Notes:
  ------
  File rowmin.dat is used to indicate those latitude lines of ncar data
  used in the interpolation process.  Rowin also indicates if the line
  is outside of the NCAR grid.  This function first checks if the line
  is outside the NCAR grid.  If the line is outside,  the function returns
  to main.  If it is inside interpolation is performed.

  Warnings
  --------
  A)  The following files MUST be generated before you use this program.
	row_index.dat
	col_index.dat
	dr.dat
	dc.dat
	binary.dat:	ncar binary file for study area

  B)  Control file information (i.e. julian_day_start, julian_day_end,
      nlines,  npixels) must be identical to that used to generate
      binary file.

***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:   D. Fraser/J. Liu       
  Last update:  2000
*****************************************************************************/


#include "ncar.h"
#include "beps.h"

//void readindex(lin,rcode)
//int lin;		/* Line number */
//int short *rcode;
//{

	//FILE *fdr, *fdc, *frow, *fcol; /* Pointers to bilinear interpolation files */
//	int bytes=4;					/* Number of bytes for a float data type */
	//int ptr, size, nitems;
   

/*	Open Binary interpolation files */
//	if ((fdr= fopen(deltarowf,"rb")) ==NULL) {
  //        printf("\n Unable to open file<%s>,  exitting program ...\n\n",
		//					deltarowf);
     //     *rcode = ERROR;
	 // return;
     //   }

//	if ((fdc=fopen(deltacolf,"rb")) == NULL) {
    //      printf("\n Unable to open file<%s>,  exitting program ...\n\n",
	//						deltacolf);
    //      *rcode = ERROR;
	//  return;
     //   }

//	if ((frow=fopen(rowf,"rb")) == NULL) {
   //       printf("\n Unable to open file<%s>,  exitting program ...\n\n",
	//						rowf);
   //       *rcode = ERROR;
	//  return;
     //   }
 
	//if ((fcol=fopen(colf,"rb")) == NULL) {
     //     printf("\n Unable to open file<%s>,  exitting program ...\n\n",
	//						colf);
      //    *rcode = ERROR;
	 // return;
      //  }
 
	
/*	*** Read a line of data from the binary interpolation files ***/

/*      Set pointer to start of output data for defined bytes */
   //     ptr=lin*NPIX_CHN + pix_offset;

/*      Seek to and Read data from row and col index files */
     //   fseek(frow,ptr,0);
     //   fseek(fcol,ptr,0);
     //   fread(&row_index[0],size=1,nitems=npixels,frow);
     //   fread(&col_index[0],size=1,nitems=npixels,fcol);

/*      Set pointer to start of line data */
      //  ptr=bytes*(lin*NPIX_CHN + pix_offset);

/*      Seek to and Read data from delta row and col files */
      //  fseek(fdr,ptr,0);
     //   fseek(fdc,ptr,0);
     //   fread(&dr[0],size=4,nitems=npixels,fdr);
     //   fread(&dc[0],size=4,nitems=npixels,fdc);

/*	   Close Input, Output Files */
      // fclose(frow);
      //  fclose(fcol);
      //  fclose(fdr);
     //   fclose(fdc);
//}


