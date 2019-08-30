
#ifndef   MY_H_FILE      
#define   MY_H_FILE       
//#include "gpubeps.cuh"
#include "beps.h"

#endif 
  
void coef_c(long pix, int short lc,float coef[],struct hy_c2 HY2[])
 {

	 float lignion_leaf,lignion_fr,lignion_wd;
	 float silt1,clay1, clay_silt1;
      float ad_factor;
	 

	  coef[35]=250.0;     //Cd  库碳氮比
	  coef[36]=34.0;      //ssd 库碳氮比
	  coef[37]=34.0;      //smd 库碳氮比
	  coef[38]=34.0;      //fsd 库碳氮比
	  coef[39]=34.0;      //fmd 库碳氮比
	  coef[40]=34.0;//25.74;     //slow库碳氮比
	  coef[41]=34.0;      //sm  库碳氮比
	  coef[42]=34.0;      //m   库碳氮比
	  coef[43]=34.0;//12.0;      //Passive 库碳氮比

	  coef[46]=34.0;   //叶C:N
	  coef[47]=250.0;   //Setm C:N
	  coef[48]=34.0;


	  if(lc==16 ||lc==17 || lc==18) {
		  lignion_leaf=0.242/(0.15+0.018*0.6*coef[36])*0.5;

		  lignion_fr=0.242/(0.15+0.018*0.6*coef[37])*0.5;

		  lignion_wd=0.50*0.5;

		  ad_factor=2.0;


	  }
	  else {
	  lignion_leaf=0.242/(0.15+0.018*0.6*coef[36])*0.8;//208-04-07: 0.8 to 1.0

	  lignion_fr=0.242/(0.15+0.018*0.6*coef[37])*0.8;

	  lignion_wd=0.50;

	  ad_factor=1.6;  

	  }

//2018年4月18日 发现 热带 雨林生物量 偏低， 其他特别是寒带地区生物量偏大， 对 系数 进行了 修改 

	  switch(lc)
	  {	
	  case 4: case 10:           // evergreeen needle conifer
		  coef[0]=0.4;  	//The ratio of NPP allocated to stem
		  coef[1]=0.10;  	//The ratio of NPP allocated to coarse roots  //  2018-05-12: 0.15 to 0.1
		  coef[2]=0.3;  	//The ratio of NPP allocated to leaves
		  coef[3]=0.2;  	//The ratio of NPP allocated to fine roots    //  2018-05-12: 0.15 to 0.2
		  coef[4]=0.0125; 	//The turn over rate of stem pool  //2014-11-10 :0.025-0.02    //2018-04-18: 0.01 to 0.0125
		  coef[5]=0.0125;	//The turn over rate of coarse root pool  //2014-11-10 :0.025-0.02  // 0.01 to 0.0125
		  coef[6]=0.4;  	//The turn over rate of leaf pool
		  coef[7]=1.0;  	 //The turn over rate of fine root pool
		  break;

	  case 1: case 7: case 8:           // evergreeen broadleaf
		  coef[0]=0.4;  	//The ratio of NPP allocated to stem  //2014-09-18: 0.5-0.4
		  coef[1]=0.1;  	//The ratio of NPP allocated to coarse roots
		  coef[2]=0.3;  	//The ratio of NPP allocated to leaves
		  coef[3]=0.2;  	//The ratio of NPP allocated to fine roots    //2014-09-18: 0.1-0.2

		  coef[4]=0.03;  //The turn over rate of stem pool             //2014-11-10 :t0.04-t0.025      //IBIS : 0.04   //2018-04-18: 0.04 to 0.03
		  coef[5]=0.03;  //The turn over rate of coarse root pool       //2014-11-10 :t0.04-t0.025     //IBIS: 0.04   //0.04 to 0.03
		  coef[6]=0.4;  	//The turn over rate of leaf pool
		  coef[7]=1.0;  	 //The turn over rate of fine root pool
		  break;


	  case 5:          // decidous needle conifer
		  coef[0]=0.4;  	//The ratio of NPP allocated to stem
		  coef[1]=0.1;  	//The ratio of NPP allocated to coarse roots
		  coef[2]=0.3;  	//The ratio of NPP allocated to leaves
		  coef[3]=0.2;  	//The ratio of NPP allocated to fine roots
		  coef[4]=0.0125;  	//The turn over rate of stem pool              //2014-11-10 :0.025-0.02   2018-04-18: 0.01 to 0.0125
		  coef[5]=0.0125; 	//The turn over rate of coarse root pool       //2014-11-10 :0.025-0.02   2018-04-18  0.01 to 0.0125  
		  coef[6]=1.0;  	//The turn over rate of leaf pool
		  coef[7]=1.0;  	 //The turn over rate of fine root pool
		  break;

	  case 2: case 3:           // decidous broadleaf
		  coef[0]=0.4;  	//The ratio of NPP allocated to stem
		  coef[1]=0.1;  	//The ratio of NPP allocated to coarse roots
		  coef[2]=0.3;  	//The ratio of NPP allocated to leaves
		  coef[3]=0.2;  	//The ratio of NPP allocated to fine roots
		  coef[4]=0.02; 	//The turn over rate of stem pool                  //2014-11-10 :0.025-0.02           
		  coef[5]=0.02; 	//The turn over rate of coarse root pool            //2014-11-10 :0.025-0.02
		  coef[6]=1.0;  	//The turn over rate of leaf pool
		  coef[7]=1.0;  	 //The turn over rate of fine root pool
		  break;

	  case 6: case 9:  //mixed
		  coef[0]=0.4;  	//The ratio of NPP allocated to stem
		  coef[1]=0.1;  	//The ratio of NPP allocated to coarse roots
		  coef[2]=0.3;  	//The ratio of NPP allocated to leaves
		  coef[3]=0.2;  	//The ratio of NPP allocated to fine roots
		  coef[4]=0.02;  	//The turn over rate of stem pool           //2014-11-10 :0.06-0.025  //2018-04-10: 0.015
		  coef[5]=0.02;  	//The turn over rate of coarse root pool     //2014-11-10 :0.06-0.025 //2018-04-10: 0.015
		  coef[6]=0.75;  	//The turn over rate of leaf pool
		  coef[7]=1.0;  	 //The turn over rate of fine root pool
		  break;

	  case 11:case 12:case 14:           //shrubs
		  coef[0]=0.2;  	//The ratio of NPP allocated to stem
		  coef[1]=0.1;  	//The ratio of NPP allocated to coarse roots
		  coef[2]=0.5;  	//The ratio of NPP allocated to leaves
		  coef[3]=0.2;  	//The ratio of NPP allocated to fine roots
		  coef[4]=0.05;  	//The turn over rate of stem pool                  //2014-11-10 :0.025-0.02    //2018-03-20: 0.04 to 0.1
		  coef[5]=0.05;  	//The turn over rate of coarse root pool           //2014-11-10 :0.025-0.02    //2018-03-20: 0.04 to 0.1
		  coef[6]=1.0;  	//The turn over rate of leaf pool
		  coef[7]=1.0;  	 //The turn over rate of fine root pool
		  break;

	  case 13: //grasslands 
		  coef[0]=0.0;  	//The ratio of NPP allocated to stem
		  coef[1]=0.0;  	//The ratio of NPP allocated to coarse roots
		  coef[2]=0.5;  	//The ratio of NPP allocated to leaves     //2017-11-01 ：0.5变为0.8
		  coef[3]=0.5;  	//The ratio of NPP allocated to fine roots   //2017-11-01 ：0.5变为0.2
		  coef[4]=0.0001; //The turn over rate of stem pool
		  coef[5]=0.0001; //The turn over rate of coarse root pool
		  coef[6]=1.0;  	//The turn over rate of leaf pool
		  coef[7]=1.0;  	 //The turn over rate of fine root pool
		  break;


	  case 16: case 17: case 18: // crops
		  coef[0]=0.4;  	//The ratio of NPP allocated to stem
		  coef[1]=0.0;  	//The ratio of NPP allocated to coarse roots
		  coef[2]=0.4;  	//The ratio of NPP allocated to leaves     //2017-11-01 ：0.5变为0.8
		  coef[3]=0.2;  	//The ratio of NPP allocated to fine roots   //2017-11-01 ：0.5变为0.2
		  coef[4]=1.0; //The turn over rate of stem pool
		  coef[5]=0.0001; //The turn over rate of coarse root pool
		  coef[6]=1.0;  	//The turn over rate of leaf pool
		  coef[7]=1.0;  	 //The turn over rate of fine root pool
		  break;
	
	  default:
		  coef[0]=0.2;  	//The ratio of NPP allocated to stem
		  coef[1]=0.1;  	//The ratio of NPP allocated to coarse roots
		  coef[2]=0.5;  	//The ratio of NPP allocated to leaves
		  coef[3]=0.2;  	//The ratio of NPP allocated to fine roots
		  coef[4]=0.0125;  	//The turn over rate of stem pool          //2017-10-30：0.025 变为0.03
		  coef[5]=0.0125;  	//The turn over rate of coarse root pool   //2017-10-30：0.025 变为0.03
		  coef[6]=0.75;  	//The turn over rate of leaf pool
		  coef[7]=1.0;  	 //The turn over rate of fine root pool
		  break;
	  }


	 silt1= HY2[pix].silt;

	 clay1= HY2[pix].clay;

	 clay_silt1=(float) (clay1+silt1); 

	 coef[8]=0.010685*(exp(-3*lignion_leaf))*((1-lignion_leaf)*0.6+lignion_leaf*0.3)* ad_factor;  
	 coef[9]=0.010685*(exp(-3*lignion_leaf))*(1-lignion_leaf)*0.4* ad_factor;  
	 coef[10]=0.010685*(exp(-3*lignion_leaf))*lignion_leaf*0.7* ad_factor; 

	 coef[11]=0.0405479*0.6* ad_factor;  
	 coef[12]=0.0405479*0.4* ad_factor;  

	 coef[13]= 0.0131517*(exp(-3*lignion_fr))*((1-lignion_fr)*0.55+lignion_fr*0.3)* ad_factor;  
	 coef[14]= 0.0131517*(exp(-3*lignion_fr))*(1-lignion_fr)*0.45* ad_factor;  
	 coef[15]=  0.0131517*(exp(-3*lignion_fr))*lignion_fr*0.7* ad_factor;  


	 coef[16]=0.05068493*0.5* ad_factor;  
	 coef[17]=0.05068493*0.5* ad_factor;  




	 coef[18]= 0.0098630*(exp(-3*lignion_wd))*((1-lignion_wd)*0.55+lignion_wd*0.45)* ad_factor;  
	 coef[19]= 0.0098630*(exp(-3*lignion_wd))*(1-lignion_wd)*0.45* ad_factor;  
	 coef[20]=  0.0098630*(exp(-3*lignion_wd))*lignion_wd*0.55* ad_factor;    // coarse litter


	 coef[21]= 0.02*(1-0.75*clay_silt1)*(0.85-0.68*clay_silt1)* ad_factor; 
	 coef[22]= 0.02*(1-0.75*clay_silt1)*(0.003+0.032*clay1)* ad_factor;
	 coef[23]= 0.02*(1-0.75*clay_silt1)
		 *(1-(0.003+0.032*clay1)-(0.85-0.68*clay_silt1)-5.0/18.0*(0.01+0.04*(1-clay_silt1)))* ad_factor; 

	 coef[24]=0.0164384*0.6* ad_factor;
	 coef[25]=0.0164384*0.4* ad_factor;


//2014-04-23: 发现 Passive 库太大 ，修改其分解系数

	 //coef[29]=0.0045*0.5/365.0* ad_factor;
	 //coef[30]=0.0045*0.5/365.0* ad_factor;

	 coef[29]=0.0045*0.5/365.0* ad_factor;
	 coef[30]=0.0045*0.5/365.0* ad_factor;

	 coef[26]=0.25*0.55/365.0* ad_factor;
	 coef[27]=0.25*(0.003-0.009*clay1)/365.0* ad_factor;
	 if(coef[27]<0.0001) coef[27]=0.0001* ad_factor; 

	 // if((0.003-0.009*clay1)>0.00001) coef[28]=0.25*(1-0.55-(0.003-0.009*clay1))/365.0;
	 //else                            coef[28]=0.25*(1-0.55-0.00001)/365.0;

	 coef[28]=0.25/365.0* ad_factor-coef[26]-coef[27];

	 //kssd_a  = coef[8] ;       kssd_sm = coef[9];   kssd_s  = coef[10];

	 //ksmd_a  = coef[11];       ksmd_sm = coef[12];

	 //kfsd_a  = coef[13];       kfsd_m  = coef[14];  kfsd_s  = coef[15];

	 //kfmd_a  = coef[16];    	  kfmd_m  = coef[17];	

	 //kcd_a   = coef[18];       kcd_m   = coef[19];  kcd_s   = coef[20];

	 //km_a    = coef[21];       km_p    = coef[22];  km_s    = coef[23];

	 //ksm_a   = coef[24];       ksm_s   = coef[25]; 

	 //ks_a    = coef[26];      ks_p    = coef[27];   ks_m    = coef[28];

	 //kp_a    = coef[29];       kp_m    = coef[30];

    }

