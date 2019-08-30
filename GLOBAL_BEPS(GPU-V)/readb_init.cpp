/*************************************************************************
  Program:     	readb_init.c 
  -------
  Description:
  -----------
	Initializes parameters to their default values
 ***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:   J. Liu /D. Fraser
  Modified by: X.F. Feng
  Last update:  July 2003
*****************************************************************************/

#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 

void readb_init(float b[]) 
{
	b[1]= 25;				/* Specific leaf area m2 kg-1 C*/
	b[4]=0.0003/RTIMES;     /* b(4) Ppt Interception co. m lai-1 d-1 */   //Ô­Îª£º0.0005
	
	b[6]=0.001/RTIMES;		/* b(6) Snowmelt Temp coeff.  */
	b[7]=0.50/RTIMES;		/* b(7) water aborption */

	b[8]=0.12/RTIMES; 		/* b(8) Snowmelt Rad coeff.  */
	b[9]=0.5; 				/* b(9) spring min. lwp  -MPa		*/
	b[10]=3000.0;			/* b(10) radiat. cut. cond. threshold kJ m-2 d-1 */
	b[11]=0.0016;			/* b(11) max. canopy cond., h20  m s-1		*/
	b[12]=2.1;				/* b(12) lwp at stomatal closure  -MPa		*/
	b[13]=0.05; 			/* b(13) slope cc humidity reduction m s-1 ug-1 m-3 */
	b[14]=432;				/* b(14) photosyn light compensat.  kJ m-2 d-1	*/
	b[15]=9720;				/* b(15) 0.5 photosyn max. kJ m-2 d-1		*/
	b[16]=6.0e-04;			/* b(16) max. mesophyll cond., CO2  m s-1	*/
	b[17]=0.0;				/* b(17) min. temp. photosyn  C		*/
	b[18]=37.0;				/* b(18) max. temp. photosyn  C		*/
	b[19]=0.0008/RTIMES;    /* b(19) leaf resp co.  kg C-1 d-1 kg-1	*/
	b[20]=0.001/RTIMES;		/* b(20) stem resp co.   kg C-1 d-1 kg-1	*/
	b[21]=0.0015/RTIMES;    /* b(21) root resp co.   kg C-1 d-1 kg-1	*/
	b[22]=-30.0;			/* b(22) snowpack energy deficit C 		*/
	b[23]=4.0;				/* b(23) temp. effect mesophyll cond.		*/
	b[24]=0.003/RTIMES;	// b(24) fine root resp co.   kg C-1 d-1 kg-1  

	b[25]=2.1;				/* b(25) Q10: constant for exp. resp. for leaf	*/
	b[26]=0.015;			/* b(26) leaf N concentration   fraction	*/
	b[27]=0.2;				/* b(27) soil characteristic parameter		*/
	b[28]=1.5;				/* b(28) Q10: constant for exp. resp. for stem	*/
	b[29]=1.9;				/* b(29) Q10: constant for exp. resp. for root	*/

	b[30]=0.25;				/* b(30) leaf carbon allocation fraction	*/
	b[31]=0.35;				/* b(31) stem carbon allocation fraction	*/
	b[32]=0.4;				/* b(32) root carbon allocation fraction	*/
  
	b[40]=0.33;				/* b(40) leaf litter C turnover   fraction	*/
	b[41]=0.0;				/* b(41) stem litter C turnover   fraction	*/
	b[42]=0.4;				/* b(42) root litter C turnover   fraction	*/
	b[43]=0.35;				/* b(43) leaf growth resp    fraction		*/
	b[44]=0.3;				/* b(44) stem growth resp    fraction		*/
	b[45]=0.35;				/* b(45) root growth resp    fraction	*/
}	
