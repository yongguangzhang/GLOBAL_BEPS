/*************************************************************************
  Program:     rs.c
  --------
  Description:
  ----------- 
	Output stomatal resistance/conductance.
***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:   J. Liu   
		(Based on BIOME-BGC)   
  Last update:	May 1998
*****************************************************************************/
	

#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 


__device__ void gpurs(int long pix,float b[],float g[],float x[],float z[],float r[],float *rrr1,float *rrr2,float *rrr3,float TI,struct hy_c1 HY1[],struct hy_c2 HY2[], struct hy_c3 HY3[])
{  

/* declaration for temporary variables/functions */
	float ppfd;
 	float tavg;
	float tmin;
	float topt, tmax;
	float vpd,vpd_open,vpd_close;
 	//float psi,psi_open,psi_close;
	float m_tavg, m_psi, m_ppfd, m_vpd;//, m_tmin,  m_co2;
	float m_most;
 	float lai;
    
	float ppfd_coef;	
	float cs,cc,cbl;
	float i;	/* index for 0:sunlit leaves and 1:shaded leaves */
	float lai_under;
	  int short		  lc_p;
    float rr[4],mm[4];
    float m0;
    float min_psi;
    float aa1;

	/* assign variables that are used more than once */
  	tavg =      z[14];      //白天温度
	tmin =      z[5];
	vpd =       0.1*z[16];		/* in kPa */
	ppfd_coef=  0.01;
	


	//topt=20+(90-z[45])*0.08;                          //为纬度的绝对值：2017-11-08由0.05改为0.08 
	
	topt=22.5+z[45]*0.25;                          //2017-11-17日修改，Z[45]为多年平均温度 
	if(topt<22.5) topt=22.5;


	tmax=36.0+z[45]*0.25;                           //原来为全球40.0， 2017-11-17日修改，Z[45]为多年平均温度  
	if(tmax<36.0) tmax=36.0;                               



	lc_p=(short  int)z[23];

//========================================================================================================================================================================//
	//if(lc_p==16 || lc_p==17 || lc_p==18)    min_psi=0.2;
	//else                                    min_psi=0.05;

min_psi=0.01;




   	//psi =       1;
	lai =       z[10];

	//psi_open =  -1;			/* in -MPa */
	//psi_close = -8;
	//vpd_open =  7.5;		/* in mb */
	//vpd_close = 20.0;
	
	


   
/*
	switch(lc_p)
	{
	case 1:// conifer 
		lai_under=1.175*exp(-0.991*z[10]);
       	break;	
	case 2:  				// deciduous forest 
		lai_under=1.5;
       	break;
	case 5:  // mixe 
		lai_under=0.5*(1.5+1.175*exp(-0.991*z[10]));
		break;	
	default:
		lai_under=0.05;
  	}													
   */

lai_under=0.01;


 for(i=0;i<=1; i++)
 {      //光对气孔导度的影响 
        if (i==0) ppfd =0.5*4.55*z[36];   //中午阳叶辐射通量密度
        else      ppfd =0.5*4.55*z[37];   //中午阴叶辐射通量密度
    
    	m_ppfd = ppfd * ppfd_coef /(1.0 + (ppfd * ppfd_coef));

	/* use the tempcurve function to generate the daily average air 
	temperature conductance multiplier */

	if(tavg<0 || tavg>tmax)  m_tavg=0.0; //温度对气孔导度的影响
		
	else
	{
		if(tavg<topt)   m_tavg=log(tavg+1)/log(topt+1);
		else    		m_tavg=cos(3.1415926*(tavg-topt)/(2*(tmax-topt)));
	}

	/* soil-leaf water potential multiplier */
	 if(x[21]<=HY1[pix].WP)  m0=(float) 0.0005;  //表层土壤含水量低于凋萎系数，很干
     else m0=(float) ((x[21]-HY1[pix].WP)/(HY1[pix].PR-HY1[pix].WP));
     
	 
	 g[41]=1.0/100.0*m0;     //土壤表层阻抗   //2017年10月20日由1/150.0 改为1/100.0
	
//============================================================================================================================================/
if( lc_p==16 || lc_p==17 ||lc_p==18){

		//第一层
		if(x[21]<=HY1[pix].WP)  mm[1]= min_psi;  //very dry__min(1.0,2.5*(x[21]-wilting)/(field-wilting))
		else if(x[21]>=HY1[pix].WP && x[21]<=HY1[pix].FC)
		{
			  mm[1]= (x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP)*2
			-(x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP)*(x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP);

		//	mm[1]=(x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP);

		}
		else mm[1]=1-0.25*(x[21]-HY1[pix].FC)/(HY1[pix].PR-HY1[pix].FC);  //2013-06-13由0.4 改为0.2

		if(mm[1]<min_psi) mm[1]=min_psi;

		//第二层
		if(x[22]<=HY2[pix].WP)  mm[2]= min_psi;  //very dry__min(1.0,2.5*(x[21]-wilting)/(field-wilting))
		else if(x[22]>=HY2[pix].WP && x[22]<=HY2[pix].FC) {
			
			mm[2]= (x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP)*2
			-(x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP)*(x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP);

			//mm[2]=(x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP);
		}
		else mm[2]=1-0.25*(x[22]-HY2[pix].FC)/(HY2[pix].PR-HY2[pix].FC);  //2013-06-13由0.4 改为0.2

		if(mm[2]<min_psi) mm[2]=min_psi;



		//第三层
		if(x[23]<HY3[pix].WP)  mm[3]= min_psi;  //very dry__min(1.0,2.5*(x[21]-wilting)/(field-wilting))
		else if(x[23]>=HY3[pix].WP && x[23]<=HY3[pix].FC) {
			 mm[3]= (x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP)*2
			-(x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP)*(x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP);

			//mm[3]=(x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP);
		}
		else mm[3]=1-0.25*(x[23]-HY3[pix].FC)/(HY3[pix].PR-HY3[pix].FC);    //2013-06-13由0.4 改为0.2


		if(mm[3]<min_psi) mm[3]=min_psi;


		//rr[0]=r[1]*mm[1]+r[2]*mm[2]+r[3]*mm[3];


		rr[1]=r[1];//*mm[1]/rr[0];
		rr[2]=r[2];//*mm[2]/rr[0];
		rr[3]=r[3];//*mm[3]/rr[0];

		m_psi=mm[1]*rr[1]+mm[2]*rr[2]+mm[3]*rr[3];

		//m_psi=rr[1]+rr[2]+rr[3];

		*rrr1=rr[1];

		*rrr2=rr[2];

		*rrr3=rr[3];

		if(m_psi<0)    m_psi=0;
		if(m_psi>1.0) m_psi=1.0;

}


else if( lc_p==11 ||lc_p==12|| lc_p==13 ||lc_p==14){

		//第一层
		if(x[21]<=HY1[pix].WP)  mm[1]= min_psi;  //very dry__min(1.0,2.5*(x[21]-wilting)/(field-wilting))
		else if(x[21]>=HY1[pix].WP && x[21]<=HY1[pix].FC)
		{
			//  mm[1]= (x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP)*2
			//-(x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP)*(x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP);

			mm[1]=(x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP);

		}
		else mm[1]=1-0.25*(x[21]-HY1[pix].FC)/(HY1[pix].PR-HY1[pix].FC);  //2013-06-13由0.4 改为0.2

		if(mm[1]<min_psi) mm[1]=min_psi;

		//第二层
		if(x[22]<=HY2[pix].WP)  mm[2]= min_psi;  //very dry__min(1.0,2.5*(x[21]-wilting)/(field-wilting))
		else if(x[22]>=HY2[pix].WP && x[22]<=HY2[pix].FC) {
			
			//mm[2]= (x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP)*2
			//-(x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP)*(x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP);

			mm[2]=(x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP);
		}
		else mm[2]=1-0.25*(x[22]-HY2[pix].FC)/(HY2[pix].PR-HY2[pix].FC);  //2013-06-13由0.4 改为0.2

		if(mm[2]<min_psi) mm[2]=min_psi;



		//第三层
		if(x[23]<HY3[pix].WP)  mm[3]= min_psi;  //very dry__min(1.0,2.5*(x[21]-wilting)/(field-wilting))
		else if(x[23]>=HY3[pix].WP && x[23]<=HY3[pix].FC) {
			// mm[3]= (x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP)*2
			//-(x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP)*(x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP);

			mm[3]=(x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP);
		}
		else mm[3]=1-0.25*(x[23]-HY3[pix].FC)/(HY3[pix].PR-HY3[pix].FC);    //2013-06-13由0.4 改为0.2


		if(mm[3]<min_psi) mm[3]=min_psi;


		//rr[0]=r[1]*mm[1]+r[2]*mm[2]+r[3]*mm[3];


		rr[1]=r[1];//*mm[1]/rr[0];
		rr[2]=r[2];//*mm[2]/rr[0];
		rr[3]=r[3];//*mm[3]/rr[0];

		m_psi=mm[1]*rr[1]+mm[2]*rr[2]+mm[3]*rr[3];

		//m_psi=rr[1]+rr[2]+rr[3];

		*rrr1=rr[1];

		*rrr2=rr[2];

		*rrr3=rr[3];

		if(m_psi<0)    m_psi=0;
		if(m_psi>1.0) m_psi=1.0;

}
else{

 //第一层
if(x[21]<=HY1[pix].WP)  mm[1]= min_psi;  //very dry__min(1.0,2.5*(x[21]-wilting)/(field-wilting))
else if(x[21]>=HY1[pix].WP && x[21]<=HY1[pix].FC)
{
	//  mm[1]= (x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP)*2
		  //  -(x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP)*(x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP);

 mm[1]=(x[21]-HY1[pix].WP)/(HY1[pix].FC-HY1[pix].WP);

}
else mm[1]=1-0.25*(x[21]-HY1[pix].FC)/(HY1[pix].PR-HY1[pix].FC);  //2013-06-13由0.4 改为0.2
 
    if(mm[1]<min_psi) mm[1]=min_psi;

//第二层
if(x[22]<=HY2[pix].WP)  mm[2]= min_psi;  //very dry__min(1.0,2.5*(x[21]-wilting)/(field-wilting))
else if(x[22]>=HY2[pix].WP && x[22]<=HY2[pix].FC) {
	// mm[2]= (x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP)*2
		   //-(x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP)*(x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP);

	mm[2]=(x[22]-HY2[pix].WP)/(HY2[pix].FC-HY2[pix].WP);
}
     else mm[2]=1-0.25*(x[22]-HY2[pix].FC)/(HY2[pix].PR-HY2[pix].FC);  //2013-06-13由0.4 改为0.2
 
	if(mm[2]<min_psi) mm[2]=min_psi;



//第三层
if(x[23]<HY3[pix].WP)  mm[3]= min_psi;  //very dry__min(1.0,2.5*(x[21]-wilting)/(field-wilting))
else if(x[23]>=HY3[pix].WP && x[23]<=HY3[pix].FC) {
	 // mm[3]= (x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP)*2
		   //-(x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP)*(x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP);

 	mm[3]=(x[23]-HY3[pix].WP)/(HY3[pix].FC-HY3[pix].WP);
}
	else mm[3]=1-0.25*(x[23]-HY3[pix].FC)/(HY3[pix].PR-HY3[pix].FC);    //2013-06-13由0.4 改为0.2
 
	
	    if(mm[3]<min_psi) mm[3]=min_psi;

 
    rr[0]=r[1]*mm[1]+r[2]*mm[2]+r[3]*mm[3];


rr[1]=r[1]*mm[1]/rr[0];
rr[2]=r[2]*mm[2]/rr[0];
rr[3]=r[3]*mm[3]/rr[0];
 
m_psi=mm[1]*rr[1]+mm[2]*rr[2]+mm[3]*rr[3];

//m_psi=rr[1]+rr[2]+rr[3];

 *rrr1=rr[1];
 
 *rrr2=rr[2];
 
 *rrr3=rr[3];

 if(m_psi<0)    m_psi=0;
  if(m_psi>1.0) m_psi=1.0;

	 }

g[80]=m_psi;

//m_psi=__max(0, m_psi);
//m_psi=__min(1, m_psi);
 

	/* CO2 multiplier */
//	m_co2 = 1.0;

	/* freezing night minimum temperature multiplier */
	//if (tmin > 0.0)        /* no effect */
	//	m_tmin = 1.0;
	//else
	//if (tmin < -8.0)       /* full tmin effect */
	//	m_tmin = 0.0;
	//else                   /* partial reduction (0.0 to -8.0 C) */
	//	m_tmin = 1.0 + (0.125 * tmin);
	
	// VPD修正因子，vapor pressure deficit multiplier, vpd in Pa 
	/*if(vpd<0.2)                       m_vpd=1.0;
    else if  (vpd>=0.2 && vpd<3.5 )   m_vpd=(3.5-vpd)/3.3;     //2012-10-14将2.8改为3.5
    else                              m_vpd=0.001;
*/

vpd_open=0.0+z[45]*0.01;
if(vpd_open<0.0) vpd_open=0.0;


vpd_close=3.0+z[45]*0.02;
if(vpd_close<3.0) vpd_close=3.0;

if(vpd<=vpd_open) m_vpd=1.0;
  else if  (vpd>vpd_open && vpd<=vpd_close )   m_vpd=(vpd_close-vpd)/(vpd_close-vpd_open);     //2012-10-14将2.8改为3.5
   else                                        m_vpd=0.001;
 


//if(vpd<0.93)                       m_vpd=1.0;
//else if  (vpd>=0.93 && vpd<4.0 )   m_vpd=pow((4.0-vpd)/(4.0-0.93),0.75)+0.0001;     //2017-05-13
 //else                              m_vpd=0.001;



	/* apply all multipliers to the maximum stomatal conductance */

//=====================================================================================================================================================================================/
	m_most= m_tavg* m_psi *  m_vpd;//;//*m_co2 * m_tmin ;
 
	//	m_most=1;
	cs =b[11]*m_ppfd*m_most;
	
	/* leaf boundary-layer conductance */
	//cbl=0.08;

       cbl=0.08;
	/* leaf cuticular conductance */
	cc=0.00005;

	/* final stomatal conductance for sunlit and shaded leaves  */

     aa1=cbl*(cc+cs)/(cbl+cc+cs);
     
	 
	// if(aa1<cs) cs=aa1;
    if (i==0) g[20] =cbl*(cc+cs)/(cbl+cc+cs);   //阳叶                   //2013年9月11日改； 原方案会压低波动
	else
	          g[21] =cbl*(cc+cs)/(cbl+cc+cs);   //阴叶

    }
	/* canopy conductance */
         g[22]=g[20]*z[32]+g[21]*z[33];

	/* canopy conductance for big leaf model */
        
		ppfd=0.5*4.55*z[20];
    
		m_ppfd = ppfd * ppfd_coef / (1.0 + (ppfd * ppfd_coef));
		
		g[22]=z[10]*b[11]*m_ppfd*m_most;

	/* understory conductance */
		ppfd=0.5*4.55*(z[9]-z[20]);
    
		m_ppfd = ppfd * ppfd_coef / (1.0 + (ppfd * ppfd_coef));
		g[19]=lai_under*b[39]*b[11]*m_ppfd*m_most;
		
    return;
} 
