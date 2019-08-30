

#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 

//#include "beps.h" 
/****************************************************************
																
	   Function:	preInitials_J(NPP0, initValues[])								
	   Description: Calculate initial values of Cx, 
					and CN ratios for a specific NPP0							
																
*****************************************************************/;
__device__ void gpupreInitials(float  NPP, float lambd,float *lambddd, float * Tmean_m, float * NPPd_m,float * coef1, float * initialValues)

{
 	

   float Cw, Ccr, Cl, Cfr;
   float Ccd, Cfsd,Cfmd,Cssd,Csmd, Csm,Cm, Cs, Cp;
   float Fm;
   float CNssd, CNsmd, CNfsd,CNfmd,CNs, CNw, CNfr,CNl,CNcd,CNsm,CNm,CNp;
	/*The ratios of carbon to nitrogen of different carbon pools*/
	
	float fw, fcr, fl, ffr; 
	float kw_cd, kcr_cd;
	float kl_sl, kfr_fl;/* sl: surface litter;fl:fine root litter*/ 
	float kssd_a, kssd_sm, kssd_s, ksmd_a, ksmd_sm;
   	float kfsd_a, kfsd_m, kfsd_s, kfmd_a, kfmd_m;
    float kcd_s, kcd_a, kcd_m;
	float km_p, km_a,km_s;
	float ksm_a, ksm_s, ks_p,ks_a,ks_m, kp_a, kp_m; 	

	float a1,a2,a3,a4,a5,a6,a7,a8;
    float lam1,lam;
    //short interval=1; 
	// Calculate the initial values of C components	
   short year,day;

    float Cw0, Ccr0, Cl0, Cfr0, Ccd0,Cssd0,Csmd0,Cfsd0,Cfmd0, Cm0, Csm0, Cs0, Cp0;
    float Cw1, Ccr1, Cl1, Cfr1, Ccd1,Cssd1,Csmd1,Cfsd1,Cfmd1, Cm1, Csm1, Cs1, Cp1;
	
	float lam_u,lam_d,part1,part2;
     
	float dCcd,dCssd,dCsmd,dCfsd,dCfmd,dCsm,dCm,dCs, dCp;
	//float ratio;
   
    float Nav;
    float Sc, Sn, NP, Nt,u1,u2,u3;
  
    float Nleaf1,Nfroot1,Nwood1, NCcroot1, Ncd1,Nssd1,Nsmd1,Nfsd1,Nfmd1,Ns1,Nsm1,Nm1,Np1;
    float Nleaf0,Nfroot0,Nwood0, NCcroot0, Ncd0,Nssd0,Nsmd0,Nfsd0,Nfmd0,Ns0,Nsm0,Nm0,Np0;
	
	float Rnpp_leaf,Rnpp_wood,Rnpp_froot,Rnpp_Ccroot;
    float totalN,totalNup,totalNmin,A1,A2;//,A3,A4;
    
 
    float  b2;
	//float TN0,TN1;
	float Nmin;//,Nav0;
   // float  B5,B6,A6;
    float CNl_av;

	//======================================================================//
 	
  
	Nav=0; 
     
    fw     = coef1[0];     //The ratio of NPP allocated to stem
	fcr    = coef1[1];     //The ratio of NPP allocated to coarse roots
	fl     = coef1[2];     //The ratio of NPP allocated to leaves
	ffr    = coef1[3];     //The ratio of NPP allocated to fine roots
	kw_cd  = coef1[4];     //The turn over rate of stem pool
	kcr_cd = coef1[5];     //The turn over rate of coarse root pool
	kl_sl  = coef1[6];     //The turn over rate of leaf pool
	kfr_fl = coef1[7];     //The turn over rate of fine root pool

    kssd_a = coef1[8];     //surface structural litter pool
    kssd_sm =coef1[9];      
    kssd_s = coef1[10];
	
	ksmd_a = coef1[11];   //surface metabolic litter  pool
	ksmd_sm =coef1[12];
    
    kfsd_a = coef1[13];   //soil structural litter pool
    kfsd_m = coef1[14];      
    kfsd_s = coef1[15];
	
	kfmd_a = coef1[16];  //soil metabolic litter  pool
	kfmd_m = coef1[17];	

	kcd_a =  coef1[18];   //coarse detritus litter  pool
    kcd_m =  coef1[19];
    kcd_s =  coef1[20];

	km_a =   coef1[21];   //soil microbial C pool
	km_p =   coef1[22];
    km_s =   coef1[23] ;

    ksm_a=   coef1[24];   //surface microbial C pool
	ksm_s=   coef1[25];

	ks_a =   coef1[26];  //slow C  pool
	ks_p =   coef1[27];
    ks_m =   coef1[28] ;

	kp_a =   coef1[29];   //passive C pool
    kp_m =   coef1[30];

	CNcd  =  coef1[35]; 
    CNssd =  coef1[36];
	CNsmd =  coef1[37]; 
    CNfsd =  coef1[38];
	CNfmd =  coef1[39];
	CNs   =  coef1[40];
	CNsm  =  coef1[41];
	CNm   =  coef1[42];
	CNp   =  coef1[43];
     
    
	CNl   =  coef1[46];  
	CNw   =  coef1[47];

   	lam1  =lambd;                                      //for surface pools
    lam   =lambd;  // for soil pools      
  
 	
	totalNup=0;
    Fm=0.2;//0.85-0.018*0.6*CNl;
    	 
	Cw =  (fw/kw_cd)* NPP ;
	Ccr=  (fcr/kcr_cd)* NPP;
	Cl =  (fl/kl_sl) * NPP;
	Cfr=  (ffr/kfr_fl)* NPP;
 
    Cssd=(1-Fm)*kl_sl*Cl/(kssd_a+kssd_sm+kssd_s);
    Csmd=Fm*kl_sl*Cl/(ksmd_a+ksmd_sm);
    Cfsd=(1-Fm)*kfr_fl*Cfr/(kfsd_a+kfsd_m+kfsd_s);
    Cfmd=Fm*kfr_fl*Cfr/(kfmd_a+kfmd_m);
    Ccd =((kw_cd*Cw+kcr_cd*Ccr)/(kcd_a+kcd_m+kcd_s));
    Cssd= Cssd/lam1; 
    Csmd= Csmd/lam1;
    Cfsd= Cfsd/lam; 
    Cfmd= Cfmd/lam;
    Ccd = Ccd/lam; 

   Csm = (Cssd*kssd_sm+Csmd*ksmd_sm)/(ksm_a+ksm_s);

   a1=Cfsd*(kfsd_m*lam*(ks_a*lam+ks_p*lam+ks_m*lam)*(kp_a*lam+kp_m*lam)+(kp_a*lam+kp_m*lam)*ks_m*lam*kfsd_s*lam+ks_p*lam*kfsd_s*lam*kp_m*lam);
   
   a2=Cfmd*kfmd_m*lam*(kp_a*lam+kp_m*lam)*(ks_a*lam+ks_p*lam+ks_m*lam);
   
   a3=Ccd*(kcd_m*lam*(ks_a*lam+ks_p*lam+ks_m*lam)*(kp_a*lam+kp_m*lam)+(kp_a*lam+kp_m*lam)*ks_m*lam*kcd_s*lam +kp_m*lam*ks_p*lam*kcd_s*lam);
   
   a4=Csm*(ksm_s*lam1*ks_m*lam*(kp_a*lam+kp_m*lam)+ks_p*lam*kp_m*lam*ksm_s*lam1);
   
   a5=Cssd*(ks_m*lam*kssd_s*lam1*(kp_a*lam+kp_m*lam)+ks_p*lam*kp_m*lam*kssd_s*lam1);
   
   a6=(km_a*lam+km_p*lam+km_s*lam)*(ks_a*lam+ks_p*lam+ks_m*lam)*(kp_a*lam+kp_m*lam);
   
   a7=km_s*lam*ks_m*lam*(kp_a*lam+kp_m*lam);
      
   a8=kp_m*lam*(km_s*lam*ks_p*lam+km_p*lam*(ks_a*lam+ks_p*lam+ks_m*lam)); 
  
   Cm=(a1+a2+a3+a4+a5)/(a6-a7-a8);
   Cs=(Csm*ksm_s*lam1+Cssd*kssd_s*lam1+Cfsd*kfsd_s*lam+Cm*km_s*lam+Ccd*kcd_s*lam)/
	   (ks_a*lam+ks_p*lam+ks_m*lam);    
   Cp=(ks_p*Cs+km_p*Cm)/(kp_a+kp_m);


    //将植被碳库年周转速率和土壤碳库年分解速率转换为每天的值
	kw_cd   =kw_cd/365.0; 
    kl_sl   =kl_sl/365.0;  
	kcr_cd  =kcr_cd /365.0;
 	kfr_fl  =kfr_fl/365.0;  
	
	kssd_a  = kssd_a /365.0;
    kssd_sm = kssd_sm/365.0;      
    kssd_s  = kssd_s/365.0;
    
	ksmd_a  = ksmd_a/365.0;
	ksmd_sm = ksmd_sm/365.0;
    
    kfsd_a  = kfsd_a/365.0;
    kfsd_m  = kfsd_m/365.0;      
    kfsd_s  = kfsd_s/365.0;
	
	kfmd_a  = kfmd_a/365.0;
	kfmd_m  = kfmd_m/365.0;	

	kcd_a   = kcd_a/365.0;
    kcd_m   = kcd_m/365.0;
    kcd_s   = kcd_s/365.0;

    km_a    = km_a/365.0;
    km_p    = km_p/365.0;
	km_s    = km_s/365.0;

    ksm_a   = ksm_a/365.0;
	ksm_s   = ksm_s/365.0;

	ks_a    = ks_a/365.0;
	ks_p    = ks_p/365.0;
    ks_m    = ks_m/365.0;

	kp_a    = kp_a/365.0;
    kp_m    = kp_m/365.0;

Nwood0=Cw/CNw;     NCcroot0=Ccr/CNw;    Nleaf0=Cl/CNl;     Nfroot0=Cfr/CNl;
Ncd0 =Ccd/CNcd;    Nssd0=Cssd/CNssd;   Nsmd0=Csmd/CNsmd;  Nfsd0=Cfsd/CNfsd;
Nfmd0=Cfmd/CNfmd;  Ns0  =Cs/CNs;       Nsm0 =Csm/CNsm;    Nm0  =Cm/CNm;
Np0  =Cp/CNp;      

Cw1 = Cw  ;   Ccr1= Ccr ;    Cl1=  Cl;    Cfr1= Cfr ;
Ccd1 = Ccd;   Cssd1=Cssd;    Csmd1=Csmd;  Cfsd1=Cfsd;
Cfmd1=Cfmd;   Csm1 =Csm;     Cm1  =Cm;    Cs1  =Cs;
Cp1  =Cp;
  
NPPd_m[0]=0;

     	 lam_u= lam1;   //for surface pool
         lam_d= lam1;   // for soil pool     
         Nmin=Cssd1 *lam_u *(kssd_a + kssd_sm + kssd_s)/CNssd+ 
	     Cfsd1 *lam_d *(kfsd_a + kfsd_m  + kfsd_s)/CNfsd+ 
		 Csmd1*(ksmd_a+ksmd_sm) * lam_u/CNsmd+
	     Cfmd1*(kfmd_a+kfmd_m) * lam_d/CNfmd+
		 Ccd1 * lam_d * (kcd_a +kcd_m +kcd_s)/CNcd+
		 Csm1 *  lam_u*( ksm_a+  ksm_s )/CNsm+ 
     	 Cm1 *lam_d*(km_a+km_s+km_p)/CNm +
		 Cs1*lam_d*(ks_a +ks_p+ks_m)/CNs+ 
		 Cp1 * lam_d * (kp_a+kp_m)/CNp 
	
		 -(Cssd1* kssd_sm +Csmd1 *ksmd_sm)*lam_u / CNsm 
		 -((Cfsd1* kfsd_m +Cfmd1 *kfmd_m+ Ccd1 * kcd_m)*lam_d+ (kp_m * Cp1+ks_m*Cs1)*lam_d)/ CNm
		 -(lam_u*(Csm1*ksm_s+ Cssd1*kssd_s)+lam_d*(km_s * Cm1+Cfsd1* kfsd_s+Ccd1 * kcd_s)) / CNs 
		 -lam_d *(km_p * Cm1+ ks_p * Cs1) / CNp;
 
		 //估计土壤有机氮库 
	     Nav=(coef1[0]*NPP/CNw+coef1[1]*NPP/CNw+coef1[2]*NPP/CNl+coef1[3]*NPP/CNl)*1; 
		 //单位为：g N m-2
		 //Nav0=Nav;//(coef1[0]*NPP/CNw+coef1[1]*NPP/CNw+coef1[2]*NPP/CNl+coef1[3]*NPP/CNl)*1000; 
 
	     Nt=0;//; TN0=Nav; Nav0=Nav;
	

 
	     for(year=0;year<20;year++){  //=============================================================/
           A1=0;A2=0;
         totalNup=0; totalN=0; totalNmin=0;
           CNl_av=0;
          for(day=1;day<365;day++){  //==========================================day loop==================================/
		 
			  if(day==364){
	       kw_cd  = coef1[4];     //The turn over rate of stem pool
	       kcr_cd = coef1[5];     //The turn over rate of coarse root pool
	       kl_sl  = coef1[6];     //The turn over rate of leaf pool
	       kfr_fl = coef1[7];     //The turn over rate of fine root pool
			  }
			  else{
            kw_cd   =0; 
            kl_sl   =0;  
	        kcr_cd  =0;
 	        kfr_fl  =0;  
			  }

         //  Fm=0.85-0.018*0.6*CNl;
         //  if(Fm<0.1) Fm=0.1;
         
		  Fm=0.2;
			
		   lam_u= lambddd[day];   //for surface pool
           lam_d= lambddd[day];   // for soil pool     

	       Cw0 = Cw1;    Ccr0 = Ccr1;   Cl0  = Cl1;   Cfr0 = Cfr1;
           Ccd0 =Ccd1;   Cssd0=Cssd1;   Csmd0=Csmd1;  Cfsd0=Cfsd1;   Cfmd0=Cfmd1;
           Csm0 =Csm1;   Cm0  =Cm1;     Cs0  =Cs1;    Cp0  =Cp1;
  
           part1   =(kw_cd * Cw0+kcr_cd * Ccr0);//(1+lam_d* (kcd_a + kcd_m + kcd_s));
           part2   =Ccd0 * lam_d* (kcd_a + kcd_m + kcd_s);//(1+lam_d* (kcd_a + kcd_m + kcd_s));
           dCcd =part1-part2;
           Ccd1 =Ccd0+ dCcd;   // Coarse detritus from woody and coarse root;
  
           Ncd0=Ccd0/CNcd;   //前一时刻该库的氮量
           //本时刻该库的氮量
		   Ncd1=Ccd0/CNcd+(kw_cd * Cw0/CNw+kcr_cd * Ccr0/CNw)        //(1+lam_d* (kcd_a + kcd_m + kcd_s))
	                      -Ccd0/CNcd * lam_d* (kcd_a + kcd_m + kcd_s);//(1+lam_d* (kcd_a + kcd_m + kcd_s));
              
            part1   =(1 - Fm)* kl_sl*Cl0;//(1+lam_u*(kssd_a + kssd_sm + kssd_s));
            part2   =Cssd0* lam_u * (kssd_a + kssd_sm + kssd_s);//(1+lam_u*(kssd_a + kssd_sm + kssd_s));
            dCssd   =part1-part2;
            Cssd1   = Cssd0+dCssd;                 // for surface structural litter
            Nssd0   = Cssd0/CNssd;                 //前一时刻该库的氮量  
             //本时刻该库的氮量
			Nssd1   =Cssd0/CNssd + (1 - Fm)* kl_sl*Cl0/CNl///(1+lam_u*(kssd_a + kssd_sm + kssd_s))
                     -Cssd0/CNssd* lam_u * (kssd_a + kssd_sm + kssd_s);///(1+lam_u*(kssd_a + kssd_sm + kssd_s));

             part1   =Fm* kl_sl * Cl0;///(1+lam_u*(ksmd_a + ksmd_sm));
             part2   = Csmd0* lam_u * (ksmd_a + ksmd_sm);///(1+lam_u*(ksmd_a + ksmd_sm));
             dCsmd   =part1-part2;
             Csmd1   = Csmd0+dCsmd;                 // for surface metabolic litter
             Nsmd0   =Csmd0/CNsmd;    //前一时刻该库的氮量
             //本时刻该库的氮量
			 Nsmd1   =Csmd0/CNsmd+ Fm* kl_sl * Cl0/CNl//(1+lam_u*(ksmd_a + ksmd_sm))
                     -Csmd0/CNsmd* lam_u * (ksmd_a + ksmd_sm);//(1+lam_u*(ksmd_a + ksmd_sm));

               
             part1=(1 - Fm)* kfr_fl* Cfr0;//(1+lam_d*(kfsd_a + kfsd_m + kfsd_s));
             part2=Cfsd0* lam_d * (kfsd_a + kfsd_m + kfsd_s);//(1+lam_d*(kfsd_a + kfsd_m + kfsd_s));
             dCfsd =part1-part2;
             Cfsd1= Cfsd0+dCfsd;                  //for soil strutural litter pool
  
             Nfsd0=Cfsd0/CNfsd; //前一时刻该库的氮量
             //本时刻该库的氮量
			 Nfsd1=Cfsd0/CNfsd+ (1 - Fm)* kfr_fl* Cfr0/CNfr//(1+lam_d*(kfsd_a + kfsd_m + kfsd_s))
                     - Cfsd0/CNfsd* lam_d * (kfsd_a + kfsd_m + kfsd_s);//(1+lam_d*(kfsd_a + kfsd_m + kfsd_s));

  
             part1=Fm * kfr_fl * Cfr0;//(1+lam_d * (kfmd_a + kfmd_m)); 
             part2=lam_d * (kfmd_a + kfmd_m)* Cfmd0;//(1+lam_d * (kfmd_a + kfmd_m));
             dCfmd=part1-part2;
             Cfmd1= Cfmd0+dCfmd;                   // for soil metabolic pool
             Nfmd0= Cfmd0/CNfmd;  //前一时刻该库的氮量
             //本时刻该库的氮量
			 Nfmd1= Cfmd0/CNfmd+Fm * kfr_fl * Cfr0/CNfr//(1+lam_d * (kfmd_a + kfmd_m))
                    -lam_d * (kfmd_a + kfmd_m)* Cfmd0/CNfmd;//(1+lam_d * (kfmd_a + kfmd_m));
			
             part1=lam_u*(Cssd0*kssd_sm+Csmd0*ksmd_sm);//(1+lam_u*(ksm_a+ksm_s));
             part2=lam_u*Csm0*(ksm_a+ksm_s);//(1+lam_u*(ksm_a+ksm_s));
             dCsm=part1-part2;
             Csm1=Csm0+dCsm;                       // for surface microbe pool
             Nsm0=Csm0/CNsm ;     //前一时刻该库的氮量
             //本时刻该库的氮量
			 Nsm1=Csm0/CNsm + lam_u*(Cssd0/CNssd*kssd_sm+Csmd0/CNsmd*ksmd_sm)//(1+lam_u*(ksm_a+ksm_s))
                    - lam_u*Csm0/CNsm*(ksm_a+ksm_s);//(1+lam_u*(ksm_a+ksm_s));


             part1=(lam_d * (kfsd_m * Cfsd0+kfmd_m*Cfmd0 + Ccd0 * kcd_m) +lam_d*(Cs0*ks_m+Cp0 * kp_m));//(1+lam_d*(km_a +  km_s +km_p));				
             part2=Cm0 * lam_d*(km_a +  km_s +km_p);//(1+lam_d*(km_a +  km_s +km_p));
             dCm=part1-part2;
             Cm1=Cm0+dCm;                          // for soil microbe pool
             Nm0=Cm0/CNm;                          //前一时刻该库的氮量
             //本时刻该库的氮量
             Nm1=Cm0/CNm +(lam_d * (kfsd_m * Cfsd0/CNfsd+kfmd_m*Cfmd0/CNfmd+Ccd0/CNcd* kcd_m)+lam_d*(Cs0/CNs*ks_m+Cp0/CNp* kp_m))
	                     //(1+lam_d*(km_a +  km_s +km_p))
                         -Cm0/CNm * lam_d*(km_a +  km_s +km_p);//(1+lam_d*(km_a +  km_s +km_p));
     
             part1=(lam_d*(Cm0*km_s + Ccd0 * kcd_s +Cfsd0*kfsd_s )+ lam_u* (Csm0*ksm_s + Cssd0*kssd_s))
	               ;//(1+lam_d *( ks_a + ks_p+ks_m));
             part2=Cs0* lam_d *( ks_a + ks_p+ks_m);//(1+lam_d *( ks_a + ks_p+ks_m));
             dCs=part1-part2;
             Cs1=Cs0+dCs;                          // for slow carbon pool
             Ns0=Cs0/CNs;        //前一时刻该库的氮量
              //本时刻该库的氮量
			 Ns1=Cs0/CNs+(lam_d*(Cm0/CNm*km_s + Ccd0/CNcd * kcd_s +Cfsd0/CNfsd*kfsd_s )+ lam_u* (Csm0/CNsm*ksm_s+Cssd0/CNssd*kssd_s))
	           //(1+lam_d *( ks_a + ks_p+ks_m))
                -Cs0/CNs* lam_d *( ks_a + ks_p+ks_m);//(1+lam_d *( ks_a + ks_p+ks_m));
   				
             dCp =(lam_d *( km_p * Cm0 + ks_p * Cs0))//(1+lam_d * (kp_m + kp_a ))
		                 - Cp0*lam_d * (kp_m + kp_a );//(1+lam_d * (kp_m + kp_a ));
             Cp1=Cp0+dCp;                         // for passive carbon pool.
             Np0=Cp0/CNp;  // //前一时刻该库的氮量		
             // //本时刻该库的氮量
			 Np1=Cp0/CNp+lam_d *( km_p * Cm0/CNm + ks_p * Cs0/CNs)//(1+lam_d * (kp_m + kp_a )) 
                        - Cp0/CNp*lam_d * (kp_m + kp_a );//(1+lam_d * (kp_m + kp_a ));

            A1=Ncd0+Nssd0+Nsmd0+Nfsd0+Nfmd0+Nsm0+Nm0+Ns0+Np0;  
            A2=Ncd1+Nssd1+Nsmd1+Nfsd1+Nfmd1+Nsm1+Nm1+Ns1+Np1;

  
       Nmin  =(float)(lam_u * (Cssd0* kssd_a/CNssd+Csmd0* ksmd_a/CNsmd+Csm0*ksm_a/CNsm) 
                      +lam_d*(Ccd0 *kcd_a/CNcd+Cfsd0* kfsd_a/CNfsd+Cfmd0 * kfmd_a/CNfmd 
                      +Cm0*km_a/CNm+Cs0*ks_a/CNs+Cp0*kp_a/CNp));       
   
       totalNmin=totalNmin+ Nmin;  
    	  
	 	Nav=Nav+Nmin;//+Ndep/365.0;                                  //假设每年的氮沉降速率为2.0 g C N/a
        if(Nav<0) Nav=0;  
       
		Sc=(Ccd1+Cssd1+Csmd1+Cfsd1+Cfmd1+Csm1+Cm1+Cs1+Cp1);
   	 		 
		Sn=Nav;
        b2=Sn/600; 
		if(b2>1.0) b2=1.0;
	 
 	NP=120*b2*exp(-8.0/100000*Sc);
        
		u1=40.8+0.01*Tmean_m[day]-0.002*Tmean_m[day]*Tmean_m[day];
		u2=0.738-0.002*Tmean_m[day];
        u3=97.42-2.504*log(NP);
        a3=Tmean_m[day]+273.37;
        a1=exp(u1-u3/(0.00831*a3));
		a2=1+exp((u2*a3-205.9)/(0.00831*a3));
        Nt=a1/a2*14.0/1000.0;//+2.0/365.0;   //g N/m-2/d
		totalNup=totalNup+Nt;
     
		Nav=(Nav-Nt);
        if(Nav<0)  Nav=0;
	    totalN=totalN+Nt; 

totalNup=totalNup+Nt;

       Rnpp_wood=coef1[0]; Rnpp_Ccroot=coef1[1];
       Rnpp_leaf=coef1[2]; Rnpp_froot=coef1[3]; 
 
       Cw1   = Cw0 +  ((NPPd_m[day]-NPPd_m[day-1])*Rnpp_wood- Cw0*kw_cd);//(1+kw_cd); 
       //Nwood0   = Cw0/CNw;       //+Nt*Rnpp_wood/1000.0  ;///(1+kw_cd); //+
       Nwood1   = Nwood0-Cw0*kw_cd/CNw;//(1+kw_cd) ;
       
	   Ccr1 =  Ccr0+  ((NPPd_m[day]-NPPd_m[day-1])*Rnpp_Ccroot- Ccr0*kcr_cd);//(1+kcr_cd);  
       //NCcroot0  =  Ccr0/CNw;   //+Nt*Rnpp_Ccroot/1000.0; //+
       NCcroot1  =  NCcroot0- Ccr0*kcr_cd/CNw;//(1+kcr_cd);
	   
	  //a1=(NPPd_m[day]-NPPd_m[day-1])*Rnpp_wood/CNw*(CNw/125.0);
      // a2=(NPPd_m[day]-NPPd_m[day-1])*Rnpp_Ccroot/CNw*(CNw/125.0); 
       
       a1=Rnpp_wood*(CNw/coef1[47]);
       a2=Rnpp_Ccroot*(CNw/coef1[47]); 


	   Cl1 = Cl0  +  ((NPPd_m[day]-NPPd_m[day-1])*Rnpp_leaf-Cl0*kl_sl);//(1+kl_sl);
      // Nleaf0   =Cl0/CNl;        // +Nt*Rnpp_leaf/1000.0;///(1+kl_sl);    //+
       Nleaf1   =Nleaf0- Cl0*kl_sl/CNl;//(1+kl_sl);        // +Nt*Rnpp_leaf/1000.0;///(1+kl_sl);    //+

	   Cfr1     = Cfr0+  ((NPPd_m[day]-NPPd_m[day-1])*Rnpp_froot-Cfr0*kfr_fl);//(1+kfr_fl);        
      // Nfroot0  = Cfr0/CNfr;  //+Nt*Rnpp_froot/1000.0;///(1+kfr_fl);   //+
       Nfroot1  = Nfroot0- Cfr0*kfr_fl/CNfr;//(1+kfr_fl);  
      
	  // a4=(NPPd_m[day]-NPPd_m[day-1])*Rnpp_leaf/CNl*(CNl/56);
       //a5=(NPPd_m[day]-NPPd_m[day-1])*Rnpp_froot/CNl*(CNl/56);
       
        a3=Rnpp_leaf*(CNl/coef1[46]);
        a4=Rnpp_froot*(CNl/coef1[46]);
       


	   Nwood1 =Nwood1 +Nt*a1/(a1+a2+a3+a4);
	   NCcroot1=NCcroot1+Nt*a2/(a1+a2+a3+a4);
	   Nleaf1 =Nleaf1 +Nt*a3/(a1+a2+a3+a4);
	   Nfroot1=Nfroot1+Nt*a4/(a1+a2+a3+a4);
	     
	  // A3=A1+Nwood0+NCcroot0+Nleaf0+Nfroot0;  //前一个时间步长的植被和土壤氮
	  // A4=A2+Nwood1+NCcroot1+Nleaf1+Nfroot1;  //当前一个时间步长的植被和土壤氮
     //  A5=Nwood1+NCcroot1+Nleaf1+Nfroot1; 
	   //TN1=A4+Nav;
	  // TN0=A3+Nav0;
	 
	   //Nav0=Nav;
	   CNw=(Cw1+Ccr1)/( Nwood1 +NCcroot1);
       CNl=(Cl1+Cfr1)/( Nleaf1 +Nfroot1);
       CNfr=CNl;
	   
     if(CNw>(coef1[47]*2)) CNw=coef1[47]*2;
     if(CNw<(coef1[47]*0.5)) CNw=coef1[47]*0.5;
    
	 if(CNl>(coef1[46]*2)) CNl=coef1[46]*2;
     if(CNl<(coef1[46]*0.5)) CNl=coef1[46]*0.5;


       Nwood0 =Nwood1;
	   NCcroot0=NCcroot1;
	   Nleaf0 =Nleaf1;
	   Nfroot0=Nfroot1;

       CNl_av=CNl_av+CNl;

    if(Ncd1<0.00001) Ncd1=0.00001; 
 	if(Nssd1<0.00001) Nssd1=0.00001;   
    if(Nsmd1<0.00001) Nsmd1=0.00001;   

    if(Nfsd1<0.00001) Nfsd1=0.00001;   
    if(Nfmd1<0.00001) Nfmd1=0.00001;   
	if(Nsm1<0.00001) Nsm1=0.00001;   
	if(Nm1<0.00001) Nm1=0.00001;   
	if(Ns1<0.00001) Ns1=0.00001;   
	if(Np1<0.00001) Np1=0.00001;    


	 CNcd=Ccd1/Ncd1;    CNssd=Cssd1/Nssd1;  CNsmd=Csmd1/Nsmd1;  CNfsd=Cfsd1/Nfsd1;
     CNfmd=Cfmd1/Nfmd1; CNsm=Csm1/Nsm1;     CNm=Cm1/Nm1;        CNs=Cs1/Ns1;
     CNp=Cp1/Np1; 

 
  //fprintf(f1,"%f, %f,%8.5f,%8.5f,%8.5f,%8.5f, %8.5f, %8.5f,%8.5f,%8.5f,%8.5f,%8.5f,%8.5f, %8.5f,%8.5f, %8.5f,%8.5f\n",
  // a1,a2,Nav+A4,totalNup,totalNmin*1000,CNw,CNl ,a4,a5, CNssd,CNsmd,CNfsd, CNfmd,CNsm,CNm, CNs,  CNp); 
 
} // =======================================the end of day loop===============================

} // the end of year loop


 

	initialValues[0] =(float) NPP;
	initialValues[1] =(float) Cw1;
	initialValues[2] =(float) Ccr1;
	initialValues[3] =(float) Cl1;
	initialValues[4] =(float) Cfr1; /* 1 to 4 is for vegetation*/
	
	initialValues[5] =(float) (Ccd1);  /* cd means coarse structural detritus*/
	initialValues[6] =(float) (Cssd1); /*ssd means surface structural detritus*/
	initialValues[7] =(float) (Csmd1); /*smd means surface metabolic  detritus*/
	initialValues[8] =(float) (Cfsd1); /*fsd means fine root strutual detritus*/
	initialValues[9] =(float) (Cfmd1); /*fmd means fine root metabolic detritus*/
	
	initialValues[10] =(float) (Csm1); /* sm is surface microbe*/
	initialValues[11] =(float) (Cm1);   /* soil microbe*/
    initialValues[12] =(float) (Cs1);   /* s represents slow carbon*/ 
 	initialValues[13] =(float) (Cp1);   /* p represents passive carbon*/

	
	initialValues[14] =(float) Fm;
    initialValues[15] =(float) 0;
	initialValues[16] =(float) CNcd;
	initialValues[17] =(float) CNssd;
	initialValues[18] =(float) CNsmd;
	initialValues[19] =(float) CNfsd;
	initialValues[20] =(float) CNfmd;
	initialValues[21] =(float) CNsm;
    initialValues[22] =(float) CNm;
    initialValues[23] =(float) CNs;
    initialValues[24] =(float) CNp;

	initialValues[25] =(float) CNw;
	initialValues[26] =(float) CNl;
	initialValues[27] =(float) Nav;
	initialValues[28] =(float) CNl_av;
 
return;    
}

