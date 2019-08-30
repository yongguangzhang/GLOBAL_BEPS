/*************************************************************************
  Program:     soilresp.c
  --------
  Description:
  ----------- 
	Output soil respiration.
***************************************************************************
  CCRS (EMS/Applications Division)
  Written by:   J. Liu      
  Last update:	May 1998

*****************************************************************************/	
#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 
/* void soilresp(coef,lambda,Cstem1,Cleaf1,Ccroot1,Cfroot1,Ccd1,Csmd1,Cssd1,
	 Cfmd1 ,Cfsd1 ,Csm1 ,Cm1 ,Cs1 ,Cp1 ,initialValues) 

	 float coef[];
 float  lambda,Cstem1,Cleaf1,Ccroot1,Cfroot1,Ccd1,Csmd1,Cssd1,Cfmd1 ,Cfsd1,Csm1 ,Cm1 ,Cs1 ,Cp1;
float  initialValues[];*/




__device__ void gpusoilresp(int long pix,int short lc,float coef[],long jday,float lat,float dTmean,float dNPP,float Ndep, float Ndep0,
struct carbonpool Cpoolold[], struct CNratio CNRold[], struct carbonpool Cpoolnew[], struct CNratio CNRnew[],float lambda, float lambda2) 

{  
   	float fw, fcr, fl, ffr; 
	float kw_cd, kcr_cd;  //Stem 和Croot的周转速率
	float kl_sl, kfr_fl,kl_sl1;  //Leaf 和froot的周转速率
	float km_p, ks_p;
   	float kssd_a, kssd_sm, kssd_s, ksmd_a, ksmd_sm,kfsd_a, kfsd_m, kfsd_s, kfmd_a, kfmd_m;
	float kcd_a, kcd_m;
   	float kcd_s,ksm_a,ksm_s, km_a, km_s, ks_a, ks_m,kp_a, kp_m; 	
	float dCcd,dCssd,dCsmd,dCfsd,dCfmd,dCsm,dCm,dCs, dCp;

	float Cw0, Ccr0, Cl0, Cfr0, Ccd0,Cssd0,Csmd0,Cfsd0,Cfmd0, Cm0, Csm0, Cs0, Cp0;
    float Cw1, Ccr1, Cl1, Cfr1, Ccd1,Cssd1,Csmd1,Cfsd1,Cfmd1, Cm1, Csm1, Cs1, Cp1;
	float Nav,NT;
   	float CNssd, CNsmd, CNfsd,CNfmd,CNs, CNw, CNl,CNcd,CNsm,CNm,CNp;

	float Sc, Sn, NP, Nt,u1,u2,u3;
    float Nmin;
    float Nleaf1,Nfroot1,Nwood1, Ncroot1, Ncd1,Nssd1,Nsmd1,Nfsd1,Nfmd1,Ns1,Nsm1,Nm1,Np1;
   float Nleaf0,Nfroot0,Nwood0, Ncroot0, Ncd0,Nssd0,Nsmd0,Nfsd0,Nfmd0,Ns0,Nsm0,Nm0,Np0;
	float Rnpp_leaf,Rnpp_wood,Rnpp_froot,Rnpp_croot;
   // float A1,A3;
    float part1,part2;	
	float Fm;	
	float lam_u,lam_d;
    short  num=0;
    float a1,a2,a3,a4,a5,a6,aa;
	float b1,b2,b3,b4,b5,b6,b7,Tk,Kt,Nav1;
   // float  totalN2;
     int j;


	lam_u= lambda;                      //for surface pool
    lam_d= lambda2;  // for soil pool     
  
    fw = coef[0];  	   fcr = coef[1]; 	fl = coef[2]; 	ffr = coef[3];   //NPP分配给生物量库的比例系数  
	  a1=1-coef[4]; a2=1-coef[5]; a3=1-coef[6]; a4=1-coef[7];
     if(a1<0.0001) a1=0.0001;
	 if(a2<0.0001) a2=0.0001;
	 if(a3<0.0001) a3=0.0001;
	 if(a4<0.0001) a4=0.0001;

     



	if(lc==1 || lc==7 || lc==8 || lc==4 || lc==10){      //常绿林 
  
	a5=1.0/90.0; a6=1.0/365.0;	
	
	if(lat>0){  //==========================================北半球==================================== 
	if(jday>=60 && jday<150){   //假设叶片在245到365天掉落 
			kw_cd   =1.0-pow(a1,a6);    
			kcr_cd  =1.0-pow(a2,a6);

			kl_sl1  =1.0-pow(a3,a5);  //	kw_cd =coef[4]*1.0/365.0;  //=========//
			kfr_fl  =1.0-pow(a4,a5);  //fl represents fine root litter 
		}
		else {
			kw_cd   =1.0-pow(a1,a6); 
			kl_sl1  =0;   	
			kcr_cd  =1.0-pow(a2,a6);
			kfr_fl  =0;
		}
	}

	else
	{//====================================================================南半球================================
		if(jday>=240 && jday<330){
			kw_cd   =1.0-pow(a1,a6);   
			kcr_cd  =1.0-pow(a2,a6);
			kl_sl1  =1.0-pow(a3,a5);     //kw_cd =coef[4]*1.0/365.0;  //=========//
			kfr_fl  =1.0-pow(a4,a5);     //fl represents fine root litter 
		}
		else {
			kw_cd   =1.0-pow(a1,a6);
			kl_sl1   =0;   
			kcr_cd  =1.0-pow(a2,a6); 
			kfr_fl  =0;
		}

	} //南半球
	
		} 
	else{  //===================================================非常绿植被==========================================================================

	 a5=1.0/65.0; a6=1.0/365.0;	
	
		if(lat>0){  //==========================================北半球==================================== 
	
	if(jday>=300&& jday<365){   //假设叶片在245到365天掉落 
    kw_cd   =1.0-pow(a1,a6);    
	kcr_cd  =1.0-pow(a2,a6);

	kl_sl1  =1.0-pow(a3,a5);  //	kw_cd =coef[4]*1.0/365.0;  //=========//
	kfr_fl  =1.0-pow(a4,a5);  //fl represents fine root litter 
	}
   else {
   kw_cd   =1.0-pow(a1,a6); 
   kl_sl1  =0;   	
   kcr_cd  =1.0-pow(a2,a6);
   kfr_fl  =0;
	}
  
  	
	}

	else
    {//====================================================================南半球================================
	if(jday>=120 && jday<180){
     kw_cd   =1.0-pow(a1,a6);   
	 kcr_cd  =1.0-pow(a2,a6);
     kl_sl1  =1.0-pow(a3,a5);     //kw_cd =coef[4]*1.0/365.0;  //=========//
	 kfr_fl  =1.0-pow(a4,a5);     //fl represents fine root litter 
	}
   else {
   kw_cd   =1.0-pow(a1,a6);
   kl_sl1   =0;   
   kcr_cd  =1.0-pow(a2,a6); 
   kfr_fl  =0;
	}
	 
	} //南半球
	
	}

 
  
	if(lc==16 || lc==17 || lc==18) kl_sl=kl_sl1;    //假设农作物地上生物量的0.6被收获
	else                           kl_sl=kl_sl1;


	kssd_a  = coef[8] ;       kssd_sm = coef[9];   kssd_s  = coef[10];
    
	ksmd_a  = coef[11];       ksmd_sm = coef[12];
    
    kfsd_a  = coef[13];       kfsd_m  = coef[14];  kfsd_s  = coef[15];
	
	kfmd_a  = coef[16];    	  kfmd_m  = coef[17];	

	kcd_a   = coef[18];       kcd_m   = coef[19];  kcd_s   = coef[20];

    km_a    = coef[21];       km_p    = coef[22];  km_s    = coef[23];

    ksm_a   = coef[24];       ksm_s   = coef[25]; 
	
	ks_a    = coef[26];      ks_p    = coef[27];   ks_m    = coef[28];

	kp_a    = coef[29];       kp_m    = coef[30];
    
	//植被碳库 
	Cw0   = Cpoolold[pix].Cstem;  Ccr0 =Cpoolold[pix].Ccroot; Cl0   =Cpoolold[pix].Cleaf; Cfr0  =Cpoolold[pix].Cfroot;
    
	//土壤碳库
	Ccd0  = Cpoolold[pix].Ccd; 
	Cssd0= Cpoolold[pix].Cssd;
	Csmd0= Cpoolold[pix].Csmd; 
	Cfsd0 = Cpoolold[pix].Cfsd; 
	Cfmd0 = Cpoolold[pix].Cfmd;
	Csm0 =  Cpoolold[pix].Csm; 
	Cm0   = Cpoolold[pix]. Cm; 
	Cs0  =  Cpoolold[pix].Cs;
	Cp0   = Cpoolold[pix].Cp;

CNcd=CNRold[pix].CNcd;
CNssd=CNRold[pix].CNssd; 
CNsmd= CNRold[pix].CNsmd;
CNfsd= CNRold[pix].CNfsd;
CNfmd= CNRold[pix].CNfmd;
CNsm = CNRold[pix].CNsm;
CNm  = CNRold[pix].CNm;
CNs  = CNRold[pix].CNs;
CNp  = CNRold[pix].CNp;

CNw  = CNRold[pix].CNstem;
CNl  = CNRold[pix].CNleaf;
Nav  =CNRold[pix].Nav;
NT  = CNRold[pix].Nt; 

if(lc==13 || lc==16 || lc==17 || lc==18)   Fm=0.6;//0.85-0.018*0.6*CNl;   //2014年8月20日由0.2 改为0.3 
else                                       Fm=0.3;


 //a1=Cw0+Ccr0+Cl0+Cfr0+Ccd0+	Cssd0+Csmd0+Cfsd0+Cfmd0+Csm0+Cm0+Cs0+Cp0;

 
           part1   =(kw_cd * Cw0+kcr_cd * Ccr0);
           part2   =Ccd0 * lam_d* (kcd_a + kcd_m + kcd_s);
           dCcd =part1-part2;
           Ccd1 =Ccd0+ dCcd;   // Coarse detritus from woody and coarse root;
  
           Ncd0=Ccd0/CNcd;   //前一时刻该库的氮量
           //本时刻该库的氮量
		   Ncd1=Ccd0/CNcd+(kw_cd * Cw0/CNw+kcr_cd * Ccr0/CNw) 
	                      -Ccd0/CNcd * lam_d* (kcd_a + kcd_m + kcd_s);
              
            part1   =(1 - Fm)* kl_sl*Cl0;
            part2   =Cssd0* lam_u * (kssd_a + kssd_sm + kssd_s);
            dCssd   =part1-part2;
            Cssd1   = Cssd0+dCssd;                 // for surface structural litter
           
			Nssd0   = Cssd0/CNssd;                 //前一时刻该库的氮量  
             //本时刻该库的氮量
			Nssd1   =Cssd0/CNssd + (1 - Fm)* kl_sl1*Cl0/CNl
                     -Cssd0/CNssd* lam_u * (kssd_a + kssd_sm + kssd_s);

             part1   =Fm* kl_sl * Cl0;
             part2   = Csmd0* lam_u * (ksmd_a + ksmd_sm);
             dCsmd   =part1-part2;
             Csmd1   = Csmd0+dCsmd;                 // for surface metabolic litter
            
			 Nsmd0   =Csmd0/CNsmd;    //前一时刻该库的氮量
             //本时刻该库的氮量
			 Nsmd1   =Csmd0/CNsmd+ Fm* kl_sl1 * Cl0/CNl
                     -Csmd0/CNsmd* lam_u * (ksmd_a + ksmd_sm);

               
             part1=(1 - Fm)* kfr_fl* Cfr0;
             part2=Cfsd0* lam_d * (kfsd_a + kfsd_m + kfsd_s);
             dCfsd =part1-part2;
             Cfsd1= Cfsd0+dCfsd;                  //for soil strutural litter pool
  
             Nfsd0=Cfsd0/CNfsd; //前一时刻该库的氮量
             //本时刻该库的氮量
			 Nfsd1=Cfsd0/CNfsd+ (1 - Fm)* kfr_fl* Cfr0/CNl
                     - Cfsd0/CNfsd* lam_d * (kfsd_a + kfsd_m + kfsd_s);

  
             part1=Fm * kfr_fl * Cfr0; 
             part2=lam_d * (kfmd_a + kfmd_m)* Cfmd0;
             dCfmd=part1-part2;
             Cfmd1= Cfmd0+dCfmd;                   // for soil metabolic pool
             Nfmd0= Cfmd0/CNfmd;  //前一时刻该库的氮量
             //本时刻该库的氮量
			 Nfmd1= Cfmd0/CNfmd+Fm * kfr_fl * Cfr0/CNl
                    -lam_d * (kfmd_a + kfmd_m)* Cfmd0/CNfmd;

             part1=lam_u*(Cssd0*kssd_sm+Csmd0*ksmd_sm);
             part2=lam_u*Csm0*(ksm_a+ksm_s);
             dCsm=part1-part2;
             Csm1=Csm0+dCsm;                       // for surface microbe pool
            
			 Nsm0=Csm0/CNsm ;     //前一时刻该库的氮量
             //本时刻该库的氮量
			 Nsm1=Csm0/CNsm + lam_u*(Cssd0/CNssd*kssd_sm+Csmd0/CNsmd*ksmd_sm)
                    - lam_u*Csm0/CNsm*(ksm_a+ksm_s);


             part1=(lam_d * (kfsd_m * Cfsd0+kfmd_m*Cfmd0 + Ccd0* kcd_m) +lam_d*(Cs0*ks_m+Cp0 * kp_m));				
             part2=Cm0 * lam_d*(km_a +  km_s +km_p);
             dCm=part1-part2;
             Cm1=Cm0+dCm;                          // for soil microbe pool
            
			 Nm0=Cm0/CNm;                          //前一时刻该库的氮量
             //本时刻该库的氮量
             Nm1=Cm0/CNm +(lam_d * (kfsd_m * Cfsd0/CNfsd+kfmd_m*Cfmd0/CNfmd+Ccd0/CNcd* kcd_m)+
				     lam_d*(Cs0/CNs*ks_m+Cp0/CNp* kp_m))
	                      -Cm0/CNm * lam_d*(km_a +  km_s +km_p);
     
         
			 part1=(lam_d*(Cm0*km_s + Ccd0 * kcd_s +Cfsd0*kfsd_s )+ lam_u* (Csm0*ksm_s + Cssd0*kssd_s));
             part2=Cs0* lam_d *( ks_a + ks_p+ks_m)/(1+lam_d *( ks_a + ks_p+ks_m));
             dCs=part1-part2;
             Cs1=Cs0+dCs;                          // for slow carbon pool
            
			 Ns0=Cs0/CNs;        //前一时刻该库的氮量
              //本时刻该库的氮量
			 Ns1=Cs0/CNs+(lam_d*(Cm0/CNm*km_s + Ccd0/CNcd * kcd_s +Cfsd0/CNfsd*kfsd_s )+ 
				 lam_u* (Csm0/CNsm*ksm_s+Cssd0/CNssd*kssd_s))
	             -Cs0/CNs* lam_d *( ks_a + ks_p+ks_m);
   				
             dCp =(lam_d *( km_p * Cm0 + ks_p * Cs0))
		                 - Cp0*lam_d * (kp_m + kp_a );
             Cp1=Cp0+dCp;                         // for passive carbon pool.
             Np0=Cp0/CNp;  // //前一时刻该库的氮量		
             // //本时刻该库的氮量
			 Np1=Cp0/CNp+lam_d *( km_p * Cm0/CNm + ks_p * Cs0/CNs) 
                        - Cp0/CNp*lam_d * (kp_m + kp_a );

			 
			// A1=Ncd0+Nssd0+Nsmd0+Nfsd0+Nfmd0+Nsm0+Nm0+Ns0+Np0;  //前一步长的总N量
  
           //  Nwood0   = Cw0/CNw;   Ncroot0  =  Ccr0/CNw;  Nleaf0   =Cl0/CNl;       Nfroot0  = Cfr0/CNl; 
            b6=(Cw0+Ccr0)/CNw+(Cl0+Cfr0)/CNl; 
            b7=(Cw0+Ccr0)/coef[47]+(Cl0+Cfr0)/coef[46]; 
            
           // A2=Ncd1+Nssd1+Nsmd1+Nfsd1+Nfmd1+Nsm1+Nm1+Ns1+Np1;  //当前步长的总N量 

         	 Nmin =(float)(lam_u * (Cssd0* kssd_a/CNssd+Csmd0* ksmd_a/CNsmd+Csm0*ksm_a/CNsm) 
                      +lam_d*(Ccd0 *kcd_a/CNcd+Cfsd0* kfsd_a/CNfsd+Cfmd0 * kfmd_a/CNfmd 
                      +Cm0*km_a/CNm+Cs0*ks_a/CNs+Cp0*kp_a/CNp));    //矿化的总N量
                  

	
  //===========================================================================================================================================//
	  if(Ndep>Ndep0)   Nav=Nav+Nmin+(Ndep-Ndep0)*0.8/365.0;   //假设氮沉的利用速率为60%
	  else             Nav=Nav+Nmin;
	  
	
		 Tk=dTmean+273.15;
		
       	Sc=Ccd1+Cssd1+Csmd1+Cfsd1+Cfmd1+Csm1+Cm1+Cs1;//+Cp1);//土壤总碳量 
   	 	Sn=Nav;                                            //土壤总氮量    
        a3=Sn/600.0; 
        
        aa=b7/b6;     
        if(aa>2.0) aa=2.0;
        if(aa<1)   aa=1.0; 


		if(a3>aa) a2=aa;//+0.25*log(a3);   //2017-10-12 居为民改， 合理性要检验， 目的是体现有N 增加， N吸收也上升
		else {

			a2=a3;
		}   

		NP=120*a2*exp(-8.0/100000*Sc);

		u1=40.8+0.01*dTmean-0.002*dTmean*dTmean;  //dTmean 为平均温度

		u2=0.738-0.002*dTmean;

		u3=97.42-2.504*log(NP);


        a1=exp(u1-u3/(0.00831*Tk)); //NT 分子项目

		a2=1+exp((u2*Tk-205.9)/(0.00831*Tk));  //NT分母项目


if(dTmean<15 && Sc>13000) Kt=(1+(15-dTmean)/30.0)*(1+(Sc-13000)/10000);
                     else Kt=1.0;

 
			 //   if(dNPP>0) Nt=a1/a2*Kt*14.0/1000.0;//植被每天吸收的N量，单位g N/m-2/d
	//  else        Nt=0;



Nt=a1/a2*Kt*14.0/1000.0;//植被每天吸收的N量，单位g N/m-2/d


		if(Nt>Nav) Nt=Nav;       //植物最多利用完全部的N
		
		Nav=Nav -Nt;     
       
		NT=NT+Nt; //吸收

		  Nav1=Nav;
	   
       Rnpp_wood=coef[0]; Rnpp_croot=coef[1]; //NPP 分配系数
       Rnpp_leaf=coef[2]; Rnpp_froot=coef[3]; 
 
       Cw1   = Cw0 +  (dNPP*Rnpp_wood- Cw0*kw_cd);///(1+kw_cd); 
       Nwood0   = Cw0/CNw;       //+Nt*Rnpp_wood/1000.0  ;///(1+kw_cd); //+
       Nwood1   = Cw0/CNw-Cw0*kw_cd/CNw;///(1+kw_cd) ;
       
	   Ccr1 =  Ccr0+  (dNPP*Rnpp_croot- Ccr0*kcr_cd);///(1+kcr_cd);  
       Ncroot0  =  Ccr0/CNw;   //+Nt*Rnpp_croot/1000.0; //+
       Ncroot1  =  Ccr0/CNw- Ccr0*kcr_cd/CNw;///(1+kcr_cd);
	   
       a1=Rnpp_wood*(CNw/coef[47])*(CNw/coef[47])*(CNw/coef[47])*(CNw/coef[47]);
       a2=Rnpp_croot*(CNw/coef[47])*(CNw/coef[47])*(CNw/coef[47])*(CNw/coef[47]); 
       
	   Cl1 = Cl0  +  (dNPP*Rnpp_leaf-Cl0*kl_sl1);///(1+kl_sl);
       Nleaf0   =Cl0/CNl;        // +Nt*Rnpp_leaf/1000.0;///(1+kl_sl);    //+
       Nleaf1   =Cl0/CNl- Cl0*kl_sl1/CNl;///(1+kl_sl);        // +Nt*Rnpp_leaf/1000.0;///(1+kl_sl);    //+

	   Cfr1     = Cfr0+  (dNPP*Rnpp_froot-Cfr0*kfr_fl);///(1+kfr_fl);        
       Nfroot0  = Cfr0/CNl;  //+Nt*Rnpp_froot/1000.0;///(1+kfr_fl);   //+
       Nfroot1  = Cfr0/CNl- Cfr0*kfr_fl/CNl;///(1+kfr_fl);  
      
	  a3=Rnpp_leaf*(CNl/coef[46])*(CNl/coef[46])*(CNl/coef[46])*(CNl/coef[46]);
      a4=Rnpp_froot*(CNl/coef[46])*(CNl/coef[46])*(CNl/coef[46])*(CNl/coef[46]);
       
	   b1=Nt*a1/(a1+a2+a3+a4);
	   b2=Nt*a2/(a1+a2+a3+a4);


	   b3=Nt*a3/(a1+a2+a3+a4);
	   b4=Nt*a4/(a1+a2+a3+a4);



	   b5=Nwood1+Ncroot1+Nleaf1+ Nfroot1;
	   
	   Nwood1 =Nwood1+b1;
	   Ncroot1=Ncroot1+b2;
	   Nleaf1 =Nleaf1+b3;
	   Nfroot1=Nfroot1+b4;
	
       //b6=Nwood1+Ncroot1+Nleaf1+ Nfroot1;


	 //  A3=A1+Nwood0+Ncroot0+Nleaf0+Nfroot0;  //前一个时间步长的植被和土壤氮
	  // A4=A2+Nwood1+Ncroot1+Nleaf1+Nfroot1;  //当前一个时间步长的植被和土壤氮
      // A5=Nwood1+Ncroot1+Nleaf1+Nfroot1; 
	   //A6=Nwood1+Ncroot1+Nleaf1+Nfroot1-(Nwood0+Ncroot0+Nleaf0+Nfroot0);
 
	if((Cw1+Ccr1)>0 &&( Nwood1 +Ncroot1)>0 )   CNw=(Cw1+Ccr1)/( Nwood1 +Ncroot1);
	else CNw  =CNRold[pix].CNstem;
	
	if((Cl1+Cfr1)>0 && ( Nleaf1 +Nfroot1)>0  )   CNl=(Cl1+Cfr1)/( Nleaf1 +Nfroot1);
	else CNl  =CNRold[pix].CNleaf;

  
if(CNw>(coef[47]*2)) {
Nav=Nav-(Cw1+Ccr1)/(coef[47]*2)+(Cw1+Ccr1)/CNw;

CNw=coef[47]*2;

}
if(Nav<0) Nav=0;


if(CNw<(coef[47]*0.5)) {
Nav=Nav-(Cw1+Ccr1)/(coef[47]*0.5)+(Cw1+Ccr1)/CNw;
CNw=coef[47]*0.5;
}

if(Nav<0) Nav=0;


if(CNl>(coef[46]*2)) {
		 CNl=coef[46]*2;
Nav=Nav-(Cl1+Cfr1)/(coef[46]*2)+(Cl1+Cfr1)/CNl;
CNl=coef[46]*2.0;

}
if(Nav<0) Nav=0;

if(CNl<(coef[46]*0.5)) { 
Nav=Nav-(Cl1+Cfr1)/(coef[46]*0.5)+(Cl1+Cfr1)/CNl;
CNl=coef[46]*0.5;

}

if(Nav<0) Nav=0;



   Nwood1  = Cw1/CNw;   Ncroot1  =  Ccr1/CNw;  Nleaf1  =Cl1/CNl;       Nfroot1 = Cfr1/CNl; 

   
   if(Ccd1<0.001)  Ccd1=0.001;
   if(Cssd1<0.001) Cssd1=0.001;
   if(Csmd1<0.001) Csmd1=0.001; 
   if(Cfsd1<0.001) Cfsd1=0.001;
   if(Cfmd1<0.001) Cfmd1=0.001;
   if(Csm1<0.001)  Csm1=0.001;
   if(Cm1<0.001)    Cm1=0.001; 
   if(Cs1<0.001)    Cs1=0.001;
   if(Cp1<0.001)    Cp1=0.001; 
   if(Cw1<0.001)    Cw1=0.001;
   if(Ccr1<0.001)   Ccr1=0.001; 
   if(Cfr1<0.001)  Cfr1=0.001;  
   if(Cl1<0.001)  Cl1=0.001;

    if(Ccd1>0)  CNcd=Ccd1/Ncd1;   
	else  CNcd=CNRold[pix].CNcd;
	 
	 if(Cssd1>0) CNssd=Cssd1/Nssd1;
	 else CNssd=CNRold[pix].CNssd;
	 
	 if(Csmd1>0) CNsmd=Csmd1/Nsmd1; 
	 else CNsmd=CNRold[pix].CNsmd;
	 
	 if(Cfsd1>0) CNfsd=Cfsd1/Nfsd1;
     else CNfsd=CNRold[pix].CNfsd;
	 
	 if(Cfmd1>0) CNfmd=Cfmd1/Nfmd1;
	 else CNfmd=CNRold[pix].CNfmd;
	 
	 if(Csm1>0) CNsm=Csm1/Nsm1;
	 else CNsm =CNRold[pix].CNsm;
	 
	 if(Cm1>0) CNm=Cm1/Nm1; 
	 else CNm  =CNRold[pix].CNm;
	 
	 if(Cs1>0) CNs=Cs1/Ns1;
     else CNs  =CNRold[pix].CNs;
	 
	 if(Cp1>0) CNp=Cp1/Np1; 
     else CNp  =CNRold[pix].CNp;


    if(CNcd>2*coef[35])   CNcd=2*coef[35]; 
    if(CNcd<0.5*coef[35]) CNcd=0.5*coef[35]; 

	if(CNssd>2*coef[36])   CNssd=2*coef[36]; 
    if(CNssd<0.5*coef[36]) CNssd=0.5*coef[36]; 

    if(CNsmd>2*coef[37])   CNsmd=2*coef[37]; 
    if(CNsmd<0.5*coef[37]) CNsmd=0.5*coef[37]; 

    if(CNfsd>2*coef[38])   CNfsd=2*coef[38]; 
    if(CNfsd<0.5*coef[38]) CNfsd=0.5*coef[38]; 

    if(CNfmd>2*coef[39])   CNfmd=2*coef[39]; 
    if(CNfmd<0.5*coef[39]) CNfmd=0.5*coef[39];

    if(CNs>2*coef[40])   CNs=2*coef[40]; 
    if(CNs<0.5*coef[40]) CNs=0.5*coef[40];

    if(CNsm>2*coef[41])   CNsm=2*coef[41]; 
    if(CNsm<0.5*coef[41]) CNsm=0.5*coef[41];
    
	if(CNm>2*coef[42])   CNm=2*coef[42]; 
    if(CNm<0.5*coef[42]) CNm=0.5*coef[42];

    if(CNp>2*coef[43])   CNp=2*coef[43]; 
    if(CNp<0.5*coef[43]) CNp=0.5*coef[43];

   // initialValues[27]=(float) (dCsmd+dCssd+dCfsd+dCfmd+dCcd+dCm+dCsm+dCs+dCp);
   
	Cpoolnew[pix].Hr=(float)(lam_u * (Cssd0* kssd_a+Csmd0* ksmd_a+Csm0*ksm_a) 
                      +lam_d*(Ccd0 *kcd_a+Cfsd0* kfsd_a+Cfmd0 * kfmd_a 
                      +Cm0*km_a+Cs0*ks_a+Cp0*kp_a))+(kl_sl1-kl_sl)*Cl0;    //heterotrophic respiration+农作物收获  
 

	 
Cpoolnew[pix].Cstem=Cw1; 
Cpoolnew[pix].Ccroot=Ccr1;
Cpoolnew[pix].Cleaf=Cl1;
Cpoolnew[pix].Cfroot=Cfr1;

Cpoolnew[pix].Ccd=Ccd1; 
Cpoolnew[pix].Cssd=Cssd1;
Cpoolnew[pix].Csmd=Csmd1; 
Cpoolnew[pix].Cfsd=Cfsd1; 
Cpoolnew[pix].Cfmd=Cfmd1;
Cpoolnew[pix].Csm=Csm1; 
Cpoolnew[pix].Cm=Cm1; 
Cpoolnew[pix].Cs=Cs1;
Cpoolnew[pix].Cp=Cp1;

CNRnew[pix].CNcd=CNcd;
CNRnew[pix].CNssd=CNssd; 
CNRnew[pix].CNsmd=CNsmd;
CNRnew[pix].CNfsd=CNfsd;
CNRnew[pix].CNfmd=CNfmd;
CNRnew[pix].CNsm=CNsm;
CNRnew[pix].CNm=CNm;
CNRnew[pix].CNs=CNs;
CNRnew[pix].CNp=CNp;

CNRnew[pix].CNstem=CNw;
CNRnew[pix].CNcroot=CNw;


CNRnew[pix].CNleaf=CNl;
CNRnew[pix].CNfroot=CNl;


CNRnew[pix].Nav=Nav;
CNRnew[pix].Nt=(float)(Cssd1+Csmd1+Csm1+Ccd1+Cfsd1+Cfmd1+Cm1+Cs1+Cp1)/1000.0; 


return;
   // return;
} 
