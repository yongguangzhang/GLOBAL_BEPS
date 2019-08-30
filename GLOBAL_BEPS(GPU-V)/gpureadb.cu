/*************************************************************************
  Program:     readb.c
  --------
  Description:
  -----------
 	Read biological parameters based on vegetation type.
***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:   J. Liu      
  Modified by: X.F. Feng
  Last update:  Sep. 2003
*****************************************************************************/

#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 


  __device__  void gpureadb(float *b,  int short lc_p)      ///zfm需要改

{

	switch(lc_p)
	{

/* reference about max. canopy cond. from BIOME-BGC and 王沙生的著作 */ 
	case 4:  case 10:          //  evergreen needleleaf, temperature zone // evergreen needleleaf zhoulei 4.1 
		   
	      b[1]=20;// 8.2;		// specific leaf area m2 kg-1 C, from BIOME-BGC  
	      b[11]=0.004;			// b(11) max. canopy cond., h20  m s-1    //关键参数    //2018-02-14 : 0.0040 to 0.0035
		 
		// b[12]=2.1;			// b(12) lwp at stomatal closure  -MPa  
		  b[12]=2.3;			// b(12) lwp at stomatal closure  -MPa , BIOME-BGC 
          b[16]=6.0e-04;		// b(16) max. mesophyll cond., CO2  m s-1      
          b[19]=0.0035/RTIMES;	// b(19) leaf resp co.  kg C-1 d-1 kg-1            /2017_11_06: 0,0025改为0.004
          b[20]=0.003/RTIMES;	// b(20) stem resp co.   kg C-1 d-1 kg-1           //2018-03-28: 0.0025 to 0.003
          b[21]=0.003/RTIMES;	// b(21) coarse root resp co.   kg C-1 d-1 kg-1    //2018-03-28: 0.0025 to 0.003
		  b[24]=0.0035/RTIMES;	// b(24) fine root resp co.   kg C-1 d-1 kg-1      //2018-03-28: 0.0025 to 0.003   2018-04-20: 0.003 to 0.004
		  
          b[25]=2.3;			// b(25) Q10: constant for exp. resp. for leaf	       
		  b[28]=2.0;			// b(28) Q10: constant for exp. resp. for stem	       
		  b[29]=2.0;			// b(29) Q10: constant for exp. resp. for roots	        
		  
		  b[39]=2.0;			// b(29) Ratio of gmax of overstory to understory	
          b[49]=45.0;            //Vcmax25  //2018-02-14: 50 to 45 
    	  b[50]=27.0;
       	break;


        case 5:                     // deciduous conifer    
		  b[1]=22.5;// 8.2;				
		  // b[1]= 10;              // specific leaf area m2 kg-1 C ,
		  // b[1]= 30;				// specific leaf area m2 kg-1 C , from BIIOME-BGC 
		 //  b[11]=0.0016;			// b(11) max. canopy cond., h20  m s-1         
          // b[11]=0.006;           // b(11) max. canopy cond., h20  m s-1 from BIIOME-BGC
		  b[11]=0.004;                                                                        // 20012_12_28,2018-02-14: 0.0040 to 0.0035
		  // b[12]=2.1;			// b(12) lwp at stomatal closure  -MPa  
		  b[12]=2.3;			// b(12) lwp at stomatal closure  -MPa , BIOME-BGC 
          b[16]=6.0e-04;		// b(16) max. mesophyll cond., CO2  m s-1      
       
		  b[19]=0.0030/RTIMES;	// b(19) leaf resp co.  kg C-1 d-1 kg-1        //2018-04-21:0.003 to 0.0035
          b[20]=0.0020/RTIMES;	// b(20) stem resp co.   kg C-1 d-1 kg-1        //2018-03-28: 0.0021 to 0.0025
          b[21]=0.0020/RTIMES;	// b(21) coarse root resp co.   kg C-1 d-1 kg-1 //2018-03-28: 0.0021 to 0.0025      
		  b[24]=0.0030/RTIMES;	// b(24) fine root resp co.   kg C-1 d-1 kg-1   //2018-03-28: 0.0033 to 0.003  2018-04-121: 0.0033 to 0.0035    
		  
		  b[25]=2.3;			// b(25) Q10: constant for exp. resp. for leaf	       
		  b[28]=2.0;			// b(28) Q10: constant for exp. resp. for stem	       
		  b[29]=2.0;			// b(29) Q10: constant for exp. resp. for roots	       
		  
		  b[39]=2.0;			// b(29) Ratio of gmax of overstory to understory	
          b[49]=45.0;             //Vcmax25  2018-02-14: 55 to 50    //2018-03-04: 50 to 45
		  b[50]=27.0;
       	break;


	case 1: case 7: case 8:                  /* evergreen broadleaf forest */ //zhoulei 4.1
           b[1]= 25.0;                       /* 1/2 (jane + biobgc) */         //2011年8月30由17.5 改为25.0                    
		 
		  b[11]=0.004;  // b(11) max. canopy cond., h20  m s-1, BIOME-BGC     //气孔导度：20171021由0.0065 改为 0.005   //20171201:0.0055 //2018-03-01:0.0050 to 0.0045
		 
          //b[12]=2.1;		   /* b(12) lwp at stomatal closure  -MPa         */
          b[12]=3.9;			// b(12) lwp at stomatal closure  -MPa , BIOME-BGC 
          b[16]=8.0e-04;	   /* b(16) max. mesophyll cond., CO2  m s-1      */
        
		  b[19]=0.0012/RTIMES;  /* b(19) leaf resp co.  kg C-1 d-1 kg-1        *///2018-04-14:0.009 改为0.0011 to 0.012    //2018-04-17: 0.012 to 0.1
          b[20]=0.0008/RTIMES;  /* b(20) stem resp co.   kg C-1 d-1 kg-1       *///20171120:0.008 改为0.0008
          b[21]=0.0008/RTIMES;  /* b(21) coarse root resp co.   kg C-1 d-1 kg-1     *///20171120:0.008 改为0.0008
		  b[24]=0.0012/RTIMES;  /* b(24) fine root resp co.   kg C-1 d-1 kg-1 */   //2018-04-14: 0.001 to 0.0012          //2018-04-17: 0.012 to 0.1          
		  
		  b[25]=2.0;			// b(25) Q10: constant for exp. resp. for leaf	     //2017-11-01：2.3改为2.0   
		  b[28]=1.9;			// b(28) Q10: constant for exp. resp. for stem	   //2017-11-01：2.1改为1.9
		  b[29]=1.9;			// b(29) Q10: constant for exp. resp. for roots	      //2017-11-01：2.13改为1.9 
		  
		  b[39]=0.5;		   /* b(29) Ratio of gmax of overstory to understory	*/
		
		  b[49]=30.0;          //Vcmax25   //2014-09-12 : 65 to 40     //2017-10-21改为37.5： //2018-02-14: 35.0 to 32.5 
		
		   b[50]=27.0;          //2017-10-26由30.0 改为25。0
		  
		  break;

	
	
	case 2: case 3:                        /* deciduous broadleaf forest */ 
		 // b[1]= 32.0;			   /* specific leaf area m2 kg-1 C, from BIOME-BGC */
		 //  b[1]= 26.6;			   /* specific leaf area m2 kg-1 C */	
		  b[1]=26.5; 
          
		  b[11]=0.005;                                //2012年12月28日由0.005改为0.006 // 20013_06_28:.0.0065     20171119:0.0065改为0.005 
		         
          b[12]=2.1;		   /* b(12) lwp at stomatal closure  -MPa         */
          // b[12]=2.2;			// b(12) lwp at stomatal closure  -MPa , BIOME-BGC 
          b[16]=8.0e-04;	   /* b(16) max. mesophyll cond., CO2  m s-1      */
        
		  b[19]=0.0035/RTIMES;  /* b(19) leaf resp co.  kg C-1 d-1 kg-1        */  //zfm5-14  原0.006 越大呼吸越大
          b[20]=0.0025/RTIMES;  /* b(20) stem resp co.   kg C-1 d-1 kg-1       *///2017-11-06: 0.0012改为0.0015
          b[21]=0.0025/RTIMES;  /* b(21) coarse root resp co.   kg C-1 d-1 kg-1       *///2017-11-06: 0.0012改为0.0015
		  b[24]=0.0035/RTIMES;  /* b(24) fine root resp co.   kg C-1 d-1 kg-1 */
		
		  b[25]=2.3;			// b(25) Q10: constant for exp. resp. for leaf	       
		  b[28]=2.0;			// b(28) Q10: constant for exp. resp. for stem	       
		  b[29]=2.0;			// b(29) Q10: constant for exp. resp. for roots	       
		  b[39]=0.5;		   /* b(29) Ratio of gmax of overstory to understory	*/
		
		  b[49]=50.0;              //2014-09-12 58 to 62               //Vcmax25  //2018-03-04: 55 to 50
		  b[50]=30.0;  
		  break;
        //////////////////////////////////////////////////////////juw
	case 6: case 9: /* mixed */ 
		  //b[1]= 17;			/* specific leaf area m2 kg-1 C */	
			b[1]= 24.0;                                      
          //b[1]= 22;			/* specific leaf area m2 kg-1 C */
		  //b[11]=0.006;		/* b(11) max. canopy cond., h20  m s-1 
		   b[11]=0.0045;		//* b(11) max. canopy cond., h20  m s-1              2018-02-14: 0.0045 to 0.0040   
           
		  b[12]=2.3;		/* b(12) lwp at stomatal closure  -MPa         */
          b[16]=7.0e-04;	/* b(16) max. mesophyll cond., CO2  m s-1      */
         
		  b[19]=0.0035/RTIMES;  // b(19) leaf resp co.  kg C-1 d-1 kg-1               2009_5_14,由0.004 改为0.0042.      2018-04-22 0.003 to 0.004
          b[20]=0.0018/RTIMES;  /* b(20) stem resp co.   kg C-1 d-1 kg-1      */      //2017-11-06: 0.0014改为0.0015     2018-04-22 0.0015 to 0.002
          b[21]=0.0018/RTIMES;  /* b(21) coarse root resp co.   kg C-1 d-1 kg-1       *///2017-11-06: 0.0014改为0.0015    2018-04-22 0.0015 to 0.002
		  b[24]=0.003/RTIMES;  /* b(24) fine root resp co.   kg C-1 d-1 kg-1       */  // 2018-04-18 0.003 to 0.0035      2018-04-22 0.003 to 0.004
		 
		  b[25]=2.3;			// b(25) Q10: constant for exp. resp. for leaf	       
		  b[28]=2.0;			// b(28) Q10: constant for exp. resp. for stem	       
		  b[29]=2.0;			// b(29) Q10: constant for exp. resp. for roots	       
		  b[39]=1.0;			/* b(29) Ratio of gmax of overstory to understory	*/
	     
		  b[49]=40.0;             //Vcmax25  2018-02-14: 45 to 40
		  b[50]=27.0;            
		  break;


/* Sep. 2003, form bio-bGC  zfm5.12*/
	case 11: case 12:   // 8:closed shrub;9:open shrub; 7: woody savanna;6:savannas; 11: permanent wetlands
        //  b[1]=12.0;				// specific leaf area m2 kg-1 C, from BIOME-BGC	
          b[1]=28.75;     //2011年7月27日由15.0 改为28.75
		  b[11]=0.0045;			                       // b(11) max. canopy cond., h20  m s-1, from BIOME-BGC  20171119:0.005改为0.0045 
		 // b[11]=0.0016;			// b(11) max. canopy cond., h20  m s-1         
		  b[12]=4.2;			// b(12) lwp at stomatal closure  -MPa , BIOME-BGC 
		 
		  b[19]=0.003/RTIMES;	// b(19) leaf resp co.  kg C-1 d-1 kg-1                                                      //2011年8月30日由0.006改为0.005
          b[20]=0.0012/RTIMES;	// b(20) stem resp co.   kg C-1 d-1 kg-1      
          b[21]=0.0012/RTIMES;	// b(21) coarse root resp co.   kg C-1 d-1 kg-1       
		  b[24]=0.003/RTIMES;	// b(24) fine root resp co.   kg C-1 d-1 kg-1                                             //2011年8月30日由0.003改为0.0025 
		 
		  b[25]=2.3;			// b(25) Q10: constant for exp. resp. for leaf	       
		  b[28]=2.1;			// b(28) Q10: constant for exp. resp. for stem	       
		  b[29]=2.1;			// b(29) Q10: constant for exp. resp. for roots	       
		  
		  b[39]=2.0;			// b(29) Ratio of gmax of overstory to understory
          b[49]=50.0;           //Vcmax25   //2018-03-04: 55 to 50
		   b[50]=27.0;
		  break; 
/* Sep. 2003, form bio-bGC  zfm5.12*/

	     case 14:   //  ============================ ================================================================ 灌木
        //  b[1]=12.0;				// specific leaf area m2 kg-1 C, from BIOME-BGC	
          b[1]=28.5;                 //2009_5_14 改为28.5； 原来15.0  ，会高估呼吸
		  b[11]=0.0045;			                           // b(11) max. canopy cond., h20  m s-1, from BIOME-BGC   20171119:0.005改为0.0035       
				  
		  // b[11]=0.0016;			// b(11) max. canopy cond., h20  m s-1         
		  b[12]=4.2;			// b(12) lwp at stomatal closure  -MPa , BIOME-BGC 
		  
		  b[19]=0.003/RTIMES;	// b(19) leaf resp co.  kg C-1 d-1 kg-1     //2011年8月30日由0.006改为0.005   
          b[20]=0.0012/RTIMES;	// b(20) stem resp co.   kg C-1 d-1 kg-1      
          b[21]=0.0012/RTIMES;	// b(21) coarse root resp co.   kg C-1 d-1 kg-1       
		  b[24]=0.003/RTIMES;	// b(24) fine root resp co.   kg C-1 d-1 kg-1    //2011年8月30日由0.003改为0.0025
		  
		  b[25]=2.3;			// b(25) Q10: constant for exp. resp. for leaf	       
		  b[28]=2.1;			// b(28) Q10: constant for exp. resp. for stem	       
		  b[29]=2.1;			// b(29) Q10: constant for exp. resp. for roots	       
		  b[39]=2.0;			// b(29) Ratio of gmax of overstory to understory 
         
		  b[49]=52.5;            //Vcmax25                               //2017-10-30：60改为50    /2017-10-30：60改为40  
		  b[50]=27.0;
		  break;

	     case 13: //  grasslands;
		   //b[1]= 10;
		   b[1]= 30;
		  //b[11]=0.019;		     // b(11) max. canopy cond., h20  m s-1  
		  b[11]=0.0045 ;              // //2018-02-14: 0.006 to 0.0045
		  b[12]=2.7;		/* b(12) lwp at stomatal closure  -MPa , from BIOME-BGC */
   		 
		  b[19]=0.004/RTIMES;	// b(19) leaf resp co.  kg C-1 d-1 kg-1     //2011年8月30日由0.006改为0.005   
		  b[20]=0.003/RTIMES;	// b(20) stem resp co.   kg C-1 d-1 kg-1      
		  b[21]=0.002/RTIMES;	// b(21) coarse root resp co.   kg C-1 d-1 kg-1       
	      b[24]=0.004/RTIMES;	// b(19) leaf resp co.  kg C-1 d-1 kg-1     //2011年8月30日由0.006改为0.005 

          b[25]=2.3;			// b(25) Q10: constant for exp. resp. for leaf	       
		  b[28]=2.1;			// b(28) Q10: constant for exp. resp. for stem	       
		  b[29]=2.1;			// b(29) Q10: constant for exp. resp. for roots	       
		  b[49]=50.0;            //Vcmax25  //2018-02-14: 55 to 50     //2018-04-15 :   50 to 55

         b[50]=30.0;
		  break;

         case 16: case 17:case 18:   // 16: cropland;
		   //b[1]= 10;
		   b[1]= 30;
		  //b[11]=0.019;		     // b(11) max. canopy cond., h20  m s-1  
		  b[11]=0.0055;                // form Jane   //5.14  原0.010 ,0.008 //5.16   //2012年12月4日由0.0060 改为0.005
		  b[12]=2.7;	             /* b(12) lwp at stomatal closure  -MPa , from BIOME-BGC */
   		
		 
		  b[19]=0.004/RTIMES;	// b(19) leaf resp co.  kg C-1 d-1 kg-1     //2011年8月30日由0.006改为0.005   
		  b[20]=0.003/RTIMES;	// b(20) stem resp co.   kg C-1 d-1 kg-1      
		  b[21]=0.0015/RTIMES;	// b(21) coarse root resp co.   kg C-1 d-1 kg-1       
	      b[24]=0.004/RTIMES;	// b(19) leaf resp co.  kg C-1 d-1 kg-1     //2011年8月30日由0.006改为0.005 

      
		  
		  b[25]=2.4;			// b(25) Q10: constant for exp. resp. for leaf	         //2017-10-29：2.3改为2.4
		  b[28]=2.1;			// b(28) Q10: constant for exp. resp. for stem	         //2017-10-29: 2.1改为2.2 
		  b[29]=2.1;			// b(29) Q10: constant for exp. resp. for roots	       

		  b[49]=60.0;                 //Vcmax25  /Vcmax25  //2012年12月28:40.0;2013_02_28:50.0          //20130623改为60.0    //2014-09-12 60 to 65  20171127:60  
		   b[50]=27.0;
		  break;

	default:       //主要是15和17类，17类的分布范围比较广


		  b[1]= 25;			    // specific leaf area m2 kg-1 C 	
       //   b[11]=0.001;		// b(11) max. canopy cond., h20  m s-1         
		  b[11]=0.0045;		    // b(11) max. canopy cond., h20  m s-1   
		   
          b[12]=2.7;		/* b(12) lwp at stomatal closure  -MPa , from BIOME-BGC */
   		  
		  b[19]=0.004/RTIMES;	// b(19) leaf resp co.  kg C-1 d-1 kg-1     //2011年8月30日由0.006改为0.005   
		  b[20]=0.0015/RTIMES;	// b(20) stem resp co.   kg C-1 d-1 kg-1      
		  b[21]=0.0015/RTIMES;	// b(21) coarse root resp co.   kg C-1 d-1 kg-1       

          b[25]=2.3;			// b(25) Q10: constant for exp. resp. for leaf	       
		  b[28]=2.1;			// b(28) Q10: constant for exp. resp. for stem	       
		  b[29]=2.1;			// b(29) Q10: constant for exp. resp. for roots	       

          b[49]=52.5;            //Vcmax25  /Vcmax25  //2012年12月28:40.0;2013_02_28:50.0       //2018-02-14: 57.5 to 52.5
		   b[50]=27.0;

  	}													
}	
