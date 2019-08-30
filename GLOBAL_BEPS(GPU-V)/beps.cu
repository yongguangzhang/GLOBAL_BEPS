
//===========================================2018-04-26-22:17


/*************************************************************************
Prograh: beps.c (Boreal Ecosystem Productivity Simulator)

===============================================================================================

BEPS model was originally developed by Chen and Liu in Canada. It is designed to calculate key
components in carbon and water cycles at regional and global scales using remotely sensed leaf
area index (LAI).This is the updated daily version of BEPS coded using GPU technique. 

The improvement of this daily version of BEPS over the earlier ones include:
(1) The carbon and nitrogen cycles are coupled. 
(2) Soil moisture dynamics is simulated using an implicit algorithm (Ju et al., 2010);

b[]  生物参数;  g[]  呼吸和蒸发散;  x[]  蒸发散;  z[]  输入要素
=================================================================================================

***********************************************************************************************/

#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 

#include "gpucoef.cu"
#include "gpudoxx.cu"
#include "gpumodel.cu"
#include "gpureadb.cu"
#include "gpupreInitial.cu"
#include "gpucarbon.cu"
#include "gpudoflux.cu"
#include "gpulw.cu"


#include "gpuinter.cu"
#include "gpumelt.cu"
#include "gpurs.cu"
#include "gpusolution.cu"
#include "gpuzcomp.cu"
#include "gpuv_moisture.cu"
#include "gpupenmon.cu"
#include "gpusoilresp.cu"
#include "gpureadb_init.cu"
#include "ncar.h"


float **dmatrix();		              // Function which allocates mem for matrix 

#define THREAD_NUM  512 //512 //512    // 512     // 512    //2017-10-02前为 512
#define BLOCK_NUM   512   // 256 //16//4         //2017-10-02前为   16        //2018-02-27: 256

#define cutilSafeCall(err) __cudaSafeCall (err, __FILE__, __LINE__)

bool InitCUDA()
{
	int count;

	cutilSafeCall(cudaGetDeviceCount(&count));
	if(count == 0) {
		fprintf(stderr, "There is no device.\n");
		return false;
	}
	int i;
	for(i = 0; i < count; i++) {
		cudaDeviceProp prop;
		cutilSafeCall(cudaGetDeviceProperties(&prop, i));
		if(cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
			if(prop.major >= 1) {
				break;
			}
		}
	}

	if(i == count) {
		fprintf(stderr, "There is no device supporting CUDA 1.x.\n");
		return false;
	}
	cutilSafeCall(cudaSetDevice(i));
	printf("%d %d\n",count,i);
	
	return true;
}

//===========================================   LAI数据转换子程序:将有效LAI转换为真实LAI ==========================

__global__ static void gpuLAI(int short*lc,float*laie,float*lai)  
{
	const int tid = threadIdx.x;
	const int bid = blockIdx.x;
 	float omega;
	int long pix;
	
	for(pix=bid * THREAD_NUM + tid;pix<npixels;pix+=BLOCK_NUM * THREAD_NUM)
	{ 
		/*
		switch(lc[pix])
		{
		case 6 : case 9://mixed =0.6, for MODIS 
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
		//==========将LAIe 转换为LAI==================================
lai[pix]=laie[pix];/// omega; //全球LAI输入数据是真实LAI,不需要再转换
	
if(lai[pix]>10.0)  lai[pix]=10.0;
if(lai[pix]<0.01)  lai[pix]=0.01;


	}          //Pixel循环结束

}

//==================计算碳库初始化需要的NPP================================================
	__global__ static void gpuModel1(long* jday, long* pix_offset,  int short*lc, float*CI,float *TI,short *Climatedata,float*lai,float*npp, xvalue*xx,						
float* lat,float*soilw1_old,float*soilw2_old,float*soilw3_old,hy_c1*HY1,hy_c2*HY2,hy_c3*HY3,float*co,float*lambd, float *lambd2)
	
	{
	float b[SIZEB],tmean,lambdt,lambdw1,lambdw2,lambdw3,x[SIZEX];
	float t_d,t_n,d_l, lambdt_d,lambdt_n;

	float soilw1_new, soilw2_new, soilw3_new;
	const int tid = threadIdx.x;
	const int bid = blockIdx.x;
	int long pix,index;
 	float CNl,coef[60];
    float a1,a2,a3;
	float w1,w2,w3;
	float lambd1_min;

   float  aa1,aa2,aa3,FC1,FC2,FC3;

	gpureadb_init(b);
		for(pix=pixoffset+bid * THREAD_NUM + tid;pix<npixels;pix+=BLOCK_NUM * THREAD_NUM)
	{ //================开始像元循环===========================

			if(lc[pix]<=0|| lc[pix]==21  || lc[pix]==19 || lai[pix]<0 ) { //对于非陆地像元，
			
			npp[pix]=0.00;  xx[pix].x11=0;  xx[pix].x21=0; xx[pix].x4=0; xx[pix].x5=0;
		}		
		else {//对于陆地像元，开始计算
			
		
			gpureadb(b,lc[pix]);	    	//Read parameters according to cover type          // read 相应的生物参数	   	
							
			index=0;
		
            gpucoef_c(pix,lc[pix],coef,HY2);
	       
		  CNl=coef[46];


			gpumodel(CNl,*jday,pix,index,lat[pix],lai[pix],lc[pix],CI[pix],TI[pix],soilw1_old[pix],soilw2_old[pix],soilw3_old[pix],
				HY1,HY2,HY3,x,b,Climatedata,xx,*co,&soilw1_new, &soilw2_new, &soilw3_new,&tmean, &t_d,&t_n,&d_l);

			npp[pix]=xx[pix].x6*1000.0;  //每天的NPP  //juw
		
			//更新3层土壤湿度=================================
			soilw1_old[pix]=soilw1_new; soilw2_old[pix]=soilw2_new; soilw3_old[pix]=soilw3_new;


			//=======计算温度对分解的调节系数， 假设土壤温度与气温相等，这是非常近似的假设
			
	

			if(t_d>35.0)           lambdt_d=1.0;
			else if( t_d<-4.0  )   lambdt_d=0.03;
			else  {
				a1=1.0/(35.0+46.32)-1/(t_d+46.32);
				a2=a1*308.56;
				lambdt_d=exp(a2);
			}

			if(t_n>35.0)        lambdt_n=1.0;
			else if( t_n<-4.0  )   lambdt_n=0.03;
			else  {
				a1=1.0/(35.0+46.32)-1/(t_n+46.32);
				a2=a1*308.56;
				lambdt_n=exp(a2);
			}


aa1=soilw1_old[pix]/HY1[pix].PR*100;   //======================================================第一层土壤实际含水量占饱和含水量的比例
FC1=HY1[pix].FC/HY1[pix].PR*100; 


aa2=soilw2_old[pix]/HY2[pix].PR*100;   //======================================================第2层土壤实际含水量占饱和含水量的比例
FC2=HY2[pix].FC/HY2[pix].PR*100; 


aa3=soilw3_old[pix]/HY3[pix].PR*100;   //======================================================第3层土壤实际含水量占饱和含水量的比例
FC3=HY3[pix].FC/HY3[pix].PR*100; 


if(FC1<10) FC1=10;    if(FC1>90) FC1=90;
if(FC2<10) FC2=10;    if(FC2>90) FC2=90;
if(FC3<10) FC3=10;    if(FC3>90) FC3=90;

gpulw1(FC1,FC2,FC3, aa1,aa2,aa3,&lambdw1,&lambdw2,&lambdw3);


/*
if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
w1=60;w2=40;w3=0;
}
else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
w1=70;w2=30;w3=0;	
}

else if(lc[pix]==5 ||lc[pix]==4){
w1=50;w2=50;w3=0;	
}

else {
w1=60;w2=40;w3=0;	
}
*/

/*
if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
	w1=70;w2=30;w3=0;
}
else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
	w1=70;w2=30;w3=0;	
}

else if(lc[pix]==5 ||lc[pix]==4){
	w1=50;w2=50;w3=0;	
}

else {
	w1=60;w2=40;w3=0;	
}
*/

if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
	w1=50;w2=40;w3=10;
}

else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
	w1=45;w2=35;w3=20;	
}

else if(lc[pix]==5 ||lc[pix]==4){
	w1=30;w2=40;w3=30;	
}

else {
	w1=35;w2=40;w3=25;	
}

a2=lambdt_d*(w1*lambdw1+w2*lambdw2+w3*lambdw3)/(w1+w2+w3);
if(a2<0.0075) a2=0.0075;

a3=lambdt_n*(w1*lambdw1+w2*lambdw2+w3*lambdw3)/(w1+w2+w3);
if(a3<0.0075) a3=0.0075;

lambd[pix]=lambd[pix]+a2*d_l/86400+a3*(1.0-d_l/86400.0);

//=========================================================================================================================================================

/*
if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
	w1=30;w2=50;w3=20;
}

else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
	w1=30;w2=50;w3=20;	


///w1=40;w2=40;w3=20;	


}

else if(lc[pix]==5 ||lc[pix]==4){
	w1=20;w2=50;w3=30;	
}

else {
	w1=40;w2=40;w3=20;	

//	w1=30;w2=50;w3=20;	

}


if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
	w1=50;w2=40;w3=10;
}

else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
	w1=55;w2=30;w3=15;	
}

else if(lc[pix]==5 ||lc[pix]==4){
	w1=30;w2=40;w3=30;	
}

else {
	w1=40;w2=40;w3=20;	
}
*/



a2=lambdt_d*(w1*lambdw1+w2*lambdw2+w3*lambdw3)/(w1+w2+w3);
if(a2<0.0075) a2=0.0075;

a3=lambdt_n*(w1*lambdw1+w2*lambdw2+w3*lambdw3)/(w1+w2+w3);
if(a3<0.0075) a3=0.0075;

lambd2[pix]=lambd2[pix]+a2*d_l/86400+a3*(1.0-d_l/86400.0);


}			   			   
	} //=========================End of PIXEL cycle PPPPPPPPPP     象元循环结束
}

//=============================================计算每天的碳通量和水通量==================================================
					
__global__ static void gpuModel2(long* jday,long* pix_offset, int short*lc,float*CI,float *TI,short * climatedata ,float *Ndep,float *Ndep0,float*lai,float*npp, xvalue*xx,float* lat,
								 float*soilw1_old,float*soilw2_old,float*soilw3_old,hy_c1*HY1,hy_c2*HY2,hy_c3*HY3,float*co,float*lambd, float *lambd2,
								 float*soilw1_n,float*soilw2_n,float*soilw3_n, carbonpool *Cpoolold,CNratio *CNRold,carbonpool *Cpoolnew,CNratio *CNRnew,
 							     float*hr,float*nep, int short *outputnep)
{
	float b[SIZEB],tmean,lambdt,lambdw1,lambdw2,lambdw3,x[SIZEX];
     float t_d,t_n,d_l,lambdt_d,lambdt_n;


	float soilw1_new, soilw2_new, soilw3_new;
	const int tid = threadIdx.x;
	const int bid = blockIdx.x;
	int long pix,index;
 	float coef[40],initialValues[40];
    int count;
    int  long  pix1; 
   float a1,a2,a3;
   float w1,w2,w3;
   float lam1,lam2, Net_dep;;

float  aa1,aa2,aa3,FC1,FC2,FC3;

	gpureadb_init(b);


	for(pix=pixoffset+bid * THREAD_NUM + tid;pix<npixels;pix+=BLOCK_NUM * THREAD_NUM)
	{ //================开始像元循环===========================

     count=pix/4950;      //4950是实际数据的列数
     pix1=   4950*count*2;

   if(lc[pix]<=0|| lc[pix]==21 || lc[pix]==19  || lai[pix]<0 ) { //对于非陆地像元，
			 
			npp[pix]=0.00; 	  xx[pix].x11=0; 	 xx[pix].x21=0; 	 xx[pix].x4=0;       xx[pix].x5=0;
		
			Cpoolnew[pix].Cstem =0.1;     Cpoolnew[pix].Cleaf =0.1;        Cpoolnew[pix].Ccroot=0.1;       Cpoolnew[pix].Cfroot=0.1; 

			Cpoolnew[pix].Ccd   = 0.1;    Cpoolnew[pix].Csmd= 0.1;         Cpoolnew[pix].Cssd= 0.1;        Cpoolnew[pix].Cfmd = 0.1;  
			
			Cpoolnew[pix].Cfsd= 0.1;      Cpoolnew[pix].Csm  = 0.1;
			Cpoolnew[pix].Cm  = 0.1;      Cpoolnew[pix].Cs   = 0.1;  Cpoolnew[pix].Cp=0.1;   
	
			
			hr[pix]    =0.0001;	nep[pix]   =0.0001;
		
            soilw1_n[pix]=0;  soilw2_n[pix]=0;soilw3_n[pix]=0;

			CNRnew[pix].CNcd   =0.1;      	CNRnew[pix].CNssd=0.1;      	CNRnew[pix].CNsmd  =0.1;

			CNRnew[pix].CNfsd  =0.1;     	CNRnew[pix].CNfmd  =0.1;     	CNRnew[pix].CNsm   =0.1;

			CNRnew[pix].CNm    =0.1;       CNRnew[pix]. CNs   =0.1;          CNRnew[pix]. CNp    =0.1;

			CNRnew[pix].CNstem =0.1;     	CNRnew[pix].CNcroot=0.1;     	CNRnew[pix].CNleaf=0.1 ;

			CNRnew[pix].CNfroot=0.1;     	CNRnew[pix].Nav     =0.1;      CNRnew[pix].  Nt     =0.1 ;


            outputnep[pix+pix1]        =-30000;      //signed short 强制把Npp变为短整型数据					
			outputnep[4950+pix+pix1]   =-30000;      //signed short 强制把Npp变为短整型数据					
			outputnep[4950*2+pix+pix1] =-30000; 
			//outputnep[4950*3+pix+pix1] =0; 
			//outputnep[4950*4+pix+pix1] =0; 
 

		}		
		else {                                            //对于陆地像元，开始计算
        
			index=1;
       	
			gpureadb(b,lc[pix]);	    	//Read parameters according to cover type          // read 相应的生物参数	   	
			
			xx[pix].x41 =Cpoolold[pix].Cleaf;   	xx[pix].x42 =Cpoolold[pix].Cstem; 	xx[pix].x43=Cpoolold[pix].Cfroot;   	xx[pix].x44=Cpoolold[pix].Ccroot;	
			

			gpumodel(CNRold[pix].CNleaf,*jday,pix,index,lat[pix],lai[pix],lc[pix],CI[pix],TI[pix],soilw1_old[pix],soilw2_old[pix],soilw3_old[pix],
			HY1,HY2,HY3,x,b, climatedata,xx,*co,&soilw1_new, &soilw2_new, &soilw3_new,&tmean,&t_d,&t_n,&d_l);

			soilw1_n[pix]=soilw1_new; soilw2_n[pix]=soilw2_new; soilw3_n[pix]=soilw3_new;

			if(lc[pix]==16 || lc[pix]==17 || lc[pix]==18){

if (soilw1_n[pix]<(HY1[pix].WP+(HY1[pix].FC-HY1[pix].WP)*0.1) )  soilw1_n[pix]=(HY1[pix].WP+(HY1[pix].FC-HY1[pix].WP)*0.1);
if (soilw2_n[pix]<(HY2[pix].WP+(HY2[pix].FC-HY2[pix].WP)*0.1) )  soilw2_n[pix]=(HY2[pix].WP+(HY2[pix].FC-HY2[pix].WP)*0.1);
if (soilw3_n[pix]<(HY3[pix].WP+(HY3[pix].FC-HY3[pix].WP)*0.1) )  soilw3_n[pix]=(HY3[pix].WP+(HY3[pix].FC-HY3[pix].WP)*0.1);

			}

			npp[pix]=xx[pix].x6*1000.0;  //每天的NPP 
			if(t_d>35.0)        lambdt_d=1.0;
			else if( t_d<-4.0  )   lambdt_d=0.03;
			else  {
				a1=1.0/(35.0+46.32)-1/(t_d+46.32);
				a2=a1*308.56;
				lambdt_d=exp(a2);
			}

			if(t_n>35.0)        lambdt_n=1.0;
			else if( t_n<-4.0  )   lambdt_n=0.03;
			else  {
				a1=1.0/(35.0+46.32)-1/(t_n+46.32);
				a2=a1*308.56;
				lambdt_n=exp(a2);
			}


			aa1=soilw1_n[pix]/HY1[pix].PR*100;   //======================================================第一层土壤实际含水量占饱和含水量的比例
			FC1=HY1[pix].FC/HY1[pix].PR*100; 


			aa2=soilw2_n[pix]/HY2[pix].PR*100;   //======================================================第2层土壤实际含水量占饱和含水量的比例
			FC2=HY2[pix].FC/HY2[pix].PR*100; 


			aa3=soilw3_n[pix]/HY3[pix].PR*100;   //======================================================第3层土壤实际含水量占饱和含水量的比例
			FC3=HY3[pix].FC/HY3[pix].PR*100; 


			if(FC1<10) FC1=10;    if(FC1>90) FC1=90;
			if(FC2<10) FC2=10;    if(FC2>90) FC2=90;
			if(FC3<10) FC3=10;    if(FC3>90) FC3=90;

			
			gpulw1(FC1,FC2,FC3, aa1,aa2,aa3,&lambdw1,&lambdw2,&lambdw3);
/*
			if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
				w1=70;w2=30;w3=0;
			}
			else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
				w1=70;w2=30;w3=0;	
			}

			else if(lc[pix]==5 ||lc[pix]==4){
				w1=50;w2=50;w3=0;	
			}

			else {
				w1=60;w2=40;w3=0;	
			}
*/


			if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
				w1=50;w2=40;w3=10;
			}

			else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
				w1=45;w2=35;w3=20;	
			}

			else if(lc[pix]==5 ||lc[pix]==4){
				w1=30;w2=40;w3=30;	
			}

			else {
				w1=35;w2=40;w3=25;	
			}

			a2=lambdt_d*(w1*lambdw1+w2*lambdw2+w3*lambdw3)/(w1+w2+w3);
			if(a2<0.0075) a2=0.0075;

			a3=lambdt_n*(w1*lambdw1+w2*lambdw2+w3*lambdw3)/(w1+w2+w3);
			if(a3<0.0075) a3=0.0075;

			lambd[pix]=a2*d_l/86400+a3*(1.0-d_l/86400.0);

			//=========================================================================================================================================================
/*
			if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
				w1=30;w2=50;w3=20;
			}

			else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
				w1=30;w2=50;w3=20;	

//w1=40;w2=40;w3=20;	


			}

			else if(lc[pix]==5 ||lc[pix]==4){
				w1=20;w2=50;w3=30;	
			}

			else {


			w1=40;w2=40;w3=20;	
			
			//w1=30;w2=50;w3=20;	
			
			}



			if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
				w1=50;w2=40;w3=10;
			}

			else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
				w1=55;w2=30;w3=15;	
			}

			else if(lc[pix]==5 ||lc[pix]==4){
				w1=30;w2=40;w3=30;	
			}

			else {
				w1=40;w2=40;w3=20;	
			}
*/

a2=lambdt_d*(w1*lambdw1+w2*lambdw2+w3*lambdw3)/(w1+w2+w3);
if(a2<0.0075) a2=0.0075;
a3=lambdt_n*(w1*lambdw1+w2*lambdw2+w3*lambdw3)/(w1+w2+w3);
if(a3<0.0075) a3=0.0075;

lambd2[pix]=a2*d_l/86400+a3*(1.0-d_l/86400.0);


//=========================================================================================================================================================
/*
if(lc[pix]==16 ||lc[pix]==17 || lc[pix]==18){
	w1=50;w2=40;w3=10;
}

else if(lc[pix]==1 ||lc[pix]==7 || lc[pix]==8){
	w1=55;w2=30;w3=15;	
}

else if(lc[pix]==5 ||lc[pix]==4){
	w1=30;w2=40;w3=30;	
}

else {
	w1=40;w2=40;w3=20;	
}
*/


			
	//==========以下代码与model1不同=================================			
	
	  gpucoef_c(pix,lc[pix],coef,HY2);



/*
     initialValues[16]=CNcd[pix] ;      	 initialValues[17]=CNssd[pix];     	 initialValues[18]=CNsmd[pix];
	
	 initialValues[19]=CNfsd[pix];      	 initialValues[20]=CNfmd[pix];   	 initialValues[21]=CNsm[pix];

	 initialValues[22]=CNm[pix];             initialValues[23]=CNs[pix];         initialValues[24]=CNp[pix];

	 initialValues[25]=CNstem[pix] ;    	 initialValues[26]=CNleaf[pix];   	 initialValues[27]=Nav[pix]  ;

	 initialValues[29]=Nt[pix]  ;
 */
	 //=================================更新碳库和不碳库的碳氮比==========================================

	 gpusoilresp(pix,lc[pix],coef,*jday,lat[pix],tmean,xx[pix].x24,Ndep[pix],Ndep0[pix],Cpoolold,CNRold,Cpoolnew,CNRnew,lambd[pix],lambd2[pix]);

/*
			Cstem_n[pix] =initialValues[1]; 			Cleaf_n[pix] =initialValues[3];          		Ccroot_n[pix]=initialValues[2];
			
			Cfroot_n[pix]=initialValues[4]; 			Ccd_n[pix]   = initialValues[5];  			Csmd_n[pix]  = initialValues[7]; 
			
			Cssd_n[pix]  = initialValues[6];			Cfmd_n[pix]  = initialValues[9];  			Cfsd_n[pix]  = initialValues[8];
			
			Csm_n[pix]   = initialValues[10];			Cm_n[pix]    = initialValues[11]; 			Cs_n[pix]    = initialValues[12];
			
			Cp_n[pix]    =initialValues[13];  

			
	CNcd_n[pix]   =initialValues[16];      	CNssd_n[pix]  =initialValues[17];      	CNsmd_n[pix]  =initialValues[18];

	CNfsd_n[pix]  =initialValues[19];     	CNfmd_n[pix]  =initialValues[20];     	CNsm_n[pix]   =initialValues[21];

    CNm_n[pix]    =initialValues[22];         CNs_n[pix]    =initialValues[23];         CNp_n[pix]    =initialValues[24];

    CNstem_n[pix] =initialValues[25];     	CNcroot_n[pix]=initialValues[25];     	CNleaf_n[pix]=initialValues[26] ;

	CNfroot_n[pix]=initialValues[26];     	Nav_n[pix]     =initialValues[27];        Nt_n[pix]     =initialValues[29] ;

*/



			hr[pix]    = hr[pix]+Cpoolnew[pix].Hr;  //异养呼吸	
			nep[pix]   = npp[pix]-hr[pix];           //NEP
			
             outputnep[pix+pix1]        =(nep[pix]*10.0); //signed short 强制把Npp变为短整型数据					
			 outputnep[4950+pix+pix1]   =(npp[pix]*10.0); //signed short 强制把Npp变为短整型数据					
			 outputnep[4950*2+pix+pix1] =(xx[pix].x11*5000.0); 
			 //outputnep[4950*3+pix+pix1] =((xx[pix].x4+xx[pix].x5)*10.0); 
			// outputnep[4950*4+pix+pix1] =((xx[pix].x31+xx[pix].x5*0)*10.0); 
}			   			   
} //=========================End of PIXEL cycle PPPPPPPPPP     象元循环结束
}

//==========================================根据gpumodel1 计算的NPP和分解速度温湿度调节因子年平均值初始碳================

__global__ static void gpuInit(long*pix_offset,hy_c2*HY2,int short*lc,float *CI,float *TI,float*lambd,float *lambd2,float*tnpp,carbonpool *Cpoolold, CNratio *CNRold)

{
	  const int tid = threadIdx.x;
	  const int bid = blockIdx.x;
 	  float coef[60];
	  int long pix;
      int long j;

      float Fm,Nav;
   
	
	 float fw, fcr, fl, ffr; 
	 float kw_cd, kcr_cd;
	 float kl_sl, kfr_fl,kl_sl1;/* sl: surface litter;fl:fine root litter*/ 
	 float kssd_a, kssd_sm, kssd_s, ksmd_a, ksmd_sm;
   	 float kfsd_a, kfsd_m, kfsd_s, kfmd_a, kfmd_m;
     float kcd_s, kcd_a, kcd_m;
	 float km_p, km_a,km_s;
	 float ksm_a, ksm_s, ks_p,ks_a,ks_m, kp_a, kp_m; 	

	float a1,a2,a3,a4,a5,a6,a7,a8;
    float lam1,lam;
    
    long year,day;

    float Cw0, Ccr0, Cl0, Cfr0, Ccd0,Cssd0,Csmd0,Cfsd0,Cfmd0, Cm0, Csm0, Cs0, Cp0;
    float Cw1, Ccr1, Cl1, Cfr1, Ccd1,Cssd1,Csmd1,Cfsd1,Cfmd1, Cm1, Csm1, Cs1, Cp1;
	
	float lam_u,lam_d,part1,part2;
     
	float dCcd,dCssd,dCsmd,dCfsd,dCfmd,dCsm,dCm,dCs, dCp;
	
    float Sc, Sn, NP, Nt;///,u1,u2,u3;
  
    float  b2;
	
    float CNl_av;

 	for(pix=pixoffset+bid * THREAD_NUM + tid;pix<npixels;pix+=BLOCK_NUM * THREAD_NUM)

	{
        if(lc[pix]>0 && lc[pix]!=21 && lc[pix]!=19) { //对于陆地像元
	
        lambd[pix]=lambd[pix]/((endyear1-startyear+1)*365.0);
          
 lambd2[pix]=lambd2[pix]/((endyear1-startyear+1)*365.0);

        gpucoef_c(pix,lc[pix],coef,HY2);
     	  
		//gpupreInitials(tnpp[pix],lambd[pix],lamb,Tmean_m,NPPd_m,coef, initValues);

  
	Nav=0; 
     
    fw     =coef[0];     //The ratio of NPP allocated to stem
	fcr    = coef[1];     //The ratio of NPP allocated to coarse roots
	fl     =coef[2];     //The ratio of NPP allocated to leaves
	ffr    = coef[3];     //The ratio of NPP allocated to fine roots
	kw_cd  = coef[4];     //The turn over rate of stem pool
	kcr_cd =coef[5];     //The turn over rate of coarse root pool
	kl_sl1  =coef[6];     //The turn over rate of leaf pool
	
	 if(lc[pix]==16 || lc[pix]==17 || lc[pix]==18) kl_sl=kl_sl1;    //假设农作物地上生物量的0.6被收获
	 else                                          kl_sl=kl_sl1;
	
	kfr_fl = coef[7];     //The turn over rate of fine root pool

    kssd_a = coef[8]*365;     //surface structural litter pool
    kssd_sm =coef[9]*365;      
    kssd_s = coef[10]*365;
	
	ksmd_a = coef[11]*365;   //surface metabolic litter  pool
	ksmd_sm =coef[12]*365;
    
    kfsd_a = coef[13]*365;   //soil structural litter pool
    kfsd_m = coef[14]*365;      
    kfsd_s = coef[15]*365;
	
	kfmd_a = coef[16]*365;  //soil metabolic litter  pool
	kfmd_m = coef[17]*365;	

	kcd_a =  coef[18]*365;   //coarse detritus litter  pool
    kcd_m =  coef[19]*365;
    kcd_s =  coef[20]*365;

	km_a =   coef[21]*365;   //soil microbial C pool
	km_p =   coef[22]*365;
    km_s =   coef[23]*365 ;

    ksm_a=   coef[24]*365;   //surface microbial C pool
	ksm_s=   coef[25]*365;

	ks_a =   coef[26]*365;  //slow C  pool
	ks_p =   coef[27]*365;
    ks_m =   coef[28] *365;

	kp_a =   coef[29]*365;   //passive C pool
    kp_m =   coef[30]*365;


   	lam1  =lambd[pix];  //for surface pools
    lam   =lambd2[pix];  // for soil pools      
  
 	
	//totalNup=0;


    if(lc[pix]==13 || lc[pix]==16 || lc[pix]==17 || lc[pix]==18) Fm=0.6;
	else 	Fm=0.3;//0.85-0.018*0.6*CNl;   //2014年8月20日由0.2 改为0.3 
    	 
	Cw1 =   (fw/kw_cd)* tnpp[pix] ;
	Ccr1=  (fcr/kcr_cd)* tnpp[pix];
	Cl1 =   (fl/kl_sl1) * tnpp[pix];
	Cfr1=  (ffr/kfr_fl)*tnpp[pix];
 
    Cssd1=(1-Fm)*kl_sl*Cl1/(kssd_a+kssd_sm+kssd_s);
    Csmd1=Fm*kl_sl*Cl1/(ksmd_a+ksmd_sm);
    Cfsd1=(1-Fm)*kfr_fl*Cfr1/(kfsd_a+kfsd_m+kfsd_s);
    Cfmd1=Fm*kfr_fl*Cfr1/(kfmd_a+kfmd_m);
    Ccd1 =((kw_cd*Cw1+kcr_cd*Ccr1)/(kcd_a+kcd_m+kcd_s));
    Cssd1= Cssd1/lam1; 
    Csmd1= Csmd1/lam1;
    Cfsd1= Cfsd1/lam; 
    Cfmd1= Cfmd1/lam;
    Ccd1 = Ccd1/lam; 

   Csm1 = (Cssd1*kssd_sm+Csmd1*ksmd_sm)/(ksm_a+ksm_s);

   a1=Cfsd1*(kfsd_m*lam*(ks_a*lam+ks_p*lam+ks_m*lam)*(kp_a*lam+kp_m*lam)+(kp_a*lam+kp_m*lam)*ks_m*lam*kfsd_s*lam+ks_p*lam*kfsd_s*lam*kp_m*lam);
   
   a2=Cfmd1*kfmd_m*lam*(kp_a*lam+kp_m*lam)*(ks_a*lam+ks_p*lam+ks_m*lam);
   
   a3=Ccd1*(kcd_m*lam*(ks_a*lam+ks_p*lam+ks_m*lam)*(kp_a*lam+kp_m*lam)+(kp_a*lam+kp_m*lam)*ks_m*lam*kcd_s*lam +kp_m*lam*ks_p*lam*kcd_s*lam);
   
   a4=Csm1*(ksm_s*lam1*ks_m*lam*(kp_a*lam+kp_m*lam)+ks_p*lam*kp_m*lam*ksm_s*lam1);
   
   a5=Cssd1*(ks_m*lam*kssd_s*lam1*(kp_a*lam+kp_m*lam)+ks_p*lam*kp_m*lam*kssd_s*lam1);
   
   a6=(km_a*lam+km_p*lam+km_s*lam)*(ks_a*lam+ks_p*lam+ks_m*lam)*(kp_a*lam+kp_m*lam);
   
   a7=km_s*lam*ks_m*lam*(kp_a*lam+kp_m*lam);
      
   a8=kp_m*lam*(km_s*lam*ks_p*lam+km_p*lam*(ks_a*lam+ks_p*lam+ks_m*lam)); 
  
   Cm1=(a1+a2+a3+a4+a5)/(a6-a7-a8);
   Cs1=(Csm1*ksm_s*lam1+Cssd1*kssd_s*lam1+Cfsd1*kfsd_s*lam+Cm1*km_s*lam+Ccd1*kcd_s*lam)/
	   (ks_a*lam+ks_p*lam+ks_m*lam);    
   Cp1=(ks_p*Cs1+km_p*Cm1)/(kp_a+kp_m);

	 
	Cpoolold[pix].Cstem =Cw1; 
	Cpoolold[pix].Cleaf =Cl1;
	Cpoolold[pix].Ccroot=Ccr1;
		
	Cpoolold[pix].Cfroot=Cfr1;
	Cpoolold[pix].Ccd   =Ccd1; 
	Cpoolold[pix].Csmd  =Csmd1 ;

	Cpoolold[pix].Cssd  =Cssd1; 
	Cpoolold[pix].Cfmd  =Cfmd1; 	
	Cpoolold[pix].Cfsd  =Cfsd1;

	Cpoolold[pix].Csm   =Csm1;
	Cpoolold[pix].Cm    =Cm1; 
	Cpoolold[pix].Cs    =Cs1;

	
	Cpoolold[pix].Cp    =Cp1;
			
    CNRold[pix].CNcd  =coef[35];     //Cd  库碳氮比
	 CNRold[pix].CNssd =coef[36];     //ssd 库碳氮比
	 CNRold[pix].CNsmd =coef[37];     //smd 库碳氮比 
	 CNRold[pix].CNfsd =coef[38];     //fsd 库碳氮比
	 CNRold[pix].CNfmd =coef[39] ;    //fmd 库碳氮比
	 CNRold[pix].CNsm  =coef[41];     //sm  库碳氮比
     CNRold[pix].CNm   =coef[42];     //sm  库碳氮比
    CNRold[pix]. CNs   =coef[40];     //slow库碳氮比
     CNRold[pix].CNp   =coef[43];     
     CNRold[pix].CNstem=coef[47];    //Setm C:NinitialValues[25] ;
	 CNRold[pix].CNcroot=coef[47];   //Setm C:NinitialValues[25] ;
	 CNRold[pix].CNleaf=coef[46] ;   //叶C:N
	 CNRold[pix].CNfroot=coef[46];   //叶C:N 
	
   Nav=(coef[0]*tnpp[pix]/CNRold[pix].CNstem+coef[1]*tnpp[pix]/CNRold[pix].CNstem+coef[2]*tnpp[pix]/CNRold[pix].CNleaf+coef[3]*tnpp[pix]/CNRold[pix].CNleaf)*20; 
	
   
   CNRold[pix].Nav     =Nav;
	
	 
	 CNRold[pix].Nt=0;	

	}

		else {

         
			Cpoolold[pix].Cstem =0.1;
			Cpoolold[pix].Cleaf =0.1;
			Cpoolold[pix].Ccroot=0.1;

			Cpoolold[pix].Cfroot=0.1;
			Cpoolold[pix].Ccd   =0.1;
			Cpoolold[pix].Csmd  =0.1;

			Cpoolold[pix].Cssd  =0.1;
			Cpoolold[pix].Cfmd  =0.1;
			Cpoolold[pix].Cfsd  =0.1;

			Cpoolold[pix].Csm   =0.1;
			Cpoolold[pix].Cm    =0.1;
			Cpoolold[pix].Cs    =0.1;


			Cpoolold[pix].Cp    =0.1;

			CNRold[pix].CNcd  =0.1;
			CNRold[pix].CNssd =0.1;
			CNRold[pix].CNsmd =0.1;
			CNRold[pix].CNfsd =0.1;
			CNRold[pix].CNfmd =0.1;
			CNRold[pix].CNsm  =0.1;
			CNRold[pix].CNm   =0.1;
			CNRold[pix]. CNs   =0.1;
			CNRold[pix].CNp   =0.1;
			CNRold[pix].CNstem=0.1;
			CNRold[pix].CNcroot=0.1;
			CNRold[pix].CNleaf=0.1;
			CNRold[pix].CNfroot=0.1;
			CNRold[pix].Nav =0.1;
			CNRold[pix].Nt=0;	

	}
	
              }

return;
}

 
//============================================================================================================================================================//
//                                                                                                                                                            //                                                                                                                                          
//                                                                                                                                                            //  
//                                                                                                                                                            //     
//                                                                                                                                                            //     
//============================================================================================================================================================//

int main()
{
    
	FILE  *ff,*ff2,*f_lat,*f_lon,*outfilenpp,*f8,*f9;


    clock_t ttt=clock();
   
    char  outfnpp_name[255],filename[200];
	unsigned char *output,*LC;
		
	int short startyearII;   
	long ptr;
    int control,control1;
	
 int short *outputnep;  
	long yearindex[200];
	int long i,j,ii;		          // zfmtry09-2-2 
	long tmpyear; 
	int short *lc;	
	float  *CI,*TI; 
	// Land cover for a line in a given study area  
	int long pix,lin;		           
	long line,pixel;                  // Position of line,pixel in China, referenced to 0,0 NOT 1,1 		  
	long jday,year,year1;             //zfm3.25
	int long lin_offset,pix_offset,factor;
    short  rcode;	    	          // Error code 
	
	short *ClimateD; 
	 
	float coef[60];                   //NPP分配比例和碳库周转速率参数数组   
	float b[SIZEB];		              // Biological parameters
	float *lai,*laie,*Ndep0,*Ndep ;	
	float *x;				          //X 变量数组  
	float *npp,*nppold, *nep,*hr;	  //NPP:net primary productivity; nep:net ecosystem productivity; hr:heterotrophic respiration
	float c[200];
	float co;
	float *tnpp,*thr;    
	
	float *lambd, *lambd2;
	     float *value;
	
  	float *soilw1_old,*soilw2_old,*soilw3_old;
    float *soilw1_n,*soilw2_n,*soilw3_n;


	float *awc;   			          // Av soil water hold cap for a line in a given study area  
	float *lon,*lat;		          // Long, lat values for a line in study area  

  
	struct xvalue *xx;		          // X variables for a line in a given study area 
	struct hy_c1 *HY1;
	struct hy_c2 *HY2;
	struct hy_c3 *HY3;
   
	//===========================================================gpu 初始化==============================================//
	//int *gpunum;
	long* gpupix_offset;
	int short *gpulc;
	float *gpuCI,*gpuTI;
	 int short *gpuoutputnep;
	short *gpusdat;
    long *gpujday;


 //   float *gpulambdd,*gpuTmean,*gpuNPPd,
	float *gpulambd,*gpulambd2,*gputnpp,*gpulai,*gpulaie,*gpunpp,*gpulat,*gpuco;
    float *gpusoilw1_old,*gpusoilw2_old,*gpusoilw3_old, *gpusoilw1_n,*gpusoilw2_n,*gpusoilw3_n;

	float *gpuNdep0,*gpuNdep;


     struct carbonpool *Cpoolold;
     struct carbonpool *Cpoolnew;

     struct CNratio  *CNRold;
     struct CNratio  *CNRnew;

float *gpuhr,*gpunep;
	
	struct xvalue *gpuxx;
	struct hy_c1 *gpuHY1;
	struct hy_c2 *gpuHY2;
	struct hy_c3 *gpuHY3;


	struct carbonpool *gpuCpoolold;
	struct carbonpool *gpuCpoolnew;

	struct CNratio  *gpuCNRold;
	struct CNratio  *gpuCNRnew;





     float ratio;
	
lin_offset=linoffset;
pix_offset=pixoffset;


	//***********============================================gpu 初始化结束======================================**********//
	 printf("CUDA initialized?\n");
	if(!InitCUDA()) {
		return 0;
	}

	printf("CUDA initialized.\n");
 	
	rcode = NOERROR;   // 0
	
	cycle_display=1;
	output_NPP='Y';
	output_GPP='Y';
	output_res='y';
	output_eva='y';
	output_tra='y';   
	output_site='n'; 
 

for(control=6;control>=1;control=control-1){   // 1: all no change; 2: Climate change only; 3 :CO2 change only; 4 : LAI change only; 5:Ndepostion 6: all change  

 //===============================================================================================================

    //-------	Read default biological paramenters  --------------------------
	readb_init(b);     //初始参数 ，45个
//========================================================================

    //读入CO2浓度数据
    ff2=fopen("\\Global_Input\\1951_2008co2.txt","rt");	
	for(year=startyear0;year<=endyear2;year++){
	fscanf(ff2,"%d %f",&ii,&c[year-startyear0]);
	printf("%d,%6.2f\n",ii, c[year-startyear0]);
    	}
    fclose(ff2);

	if((f_lat=fopen("\\Global_Input\\Global_latlon\\globallat0727.dat","rb"))==NULL) printf("Can not open latitude file\n");    
	if((f_lon=fopen("\\Global_Input\\Global_latlon\\globallon0727.dat","rb"))==NULL) printf("Can not open longitude file"); // read lat file, X.F April

	cudaMalloc((void**) &gpupix_offset, sizeof(long));     //right                                  //1
	cudaMemcpy(gpupix_offset, &pix_offset, sizeof(long), cudaMemcpyHostToDevice);

	control1=control;
	
	//==== ================================================================================================================= 开始行循环===============================  
			
	for (lin =0 ; lin <(nlines-lin_offset);lin=lin+1 ) { //===========================================================================================================
	    	
	line = lin_offset + lin;     //line 为实际行号
	//printf("line=%d  %d  %d\n",line, lin_offset,lin);

    LC=(unsigned char*)malloc(npixels*sizeof(unsigned char));   //1       //纬度, 浮点
    
	lat=(float*)malloc(npixels*sizeof(float));   //2       //纬度, 浮点
	lon=(float*)malloc(npixels*sizeof(float));      //3    //经度, 浮点
	 
	lc=(int short*)malloc(npixels*sizeof(int short));//4
	awc=(float*)malloc(npixels*sizeof(float));	  //5   
	lai=(float*)malloc(npixels*sizeof(float));  //6
	laie=(float*)malloc(npixels*sizeof(float));  //7
    
	Ndep0=(float*)malloc(npixels*sizeof(float));  //8
	Ndep=(float*)malloc(npixels*sizeof(float));  //9
	
	nppold=(float*)malloc(npixels*sizeof(float)); //10
	npp=(float*)malloc(npixels*sizeof(float));//11
	nep=(float*)malloc(npixels*sizeof(float));//12
	hr=(float*)malloc(npixels*sizeof(float));//13
	
value=(float*)malloc(npixels*sizeof(float));//13;




	output=(unsigned char*)malloc(npixels*sizeof(unsigned char));//14
	outputnep=( int short*)malloc(npixels3*sizeof( int short));//15  //X.F June 20, 2003

	tnpp = (float *)malloc(npixels*sizeof(float));//16  //X.F June 20, 2003   
   
	thr = (float *)malloc(npixels*sizeof(float));//16  //X.F June 20, 2003   


	lambd=(float *)malloc(npixels*sizeof(float));//17   //X.F June 20, 2003;  
	lambd2=(float *)malloc(npixels*sizeof(float));//17   //X.F June 20, 2003;  	
		
	Cpoolold=(struct carbonpool*)malloc(npixels*sizeof(struct carbonpool));  //49
	Cpoolnew=(struct carbonpool*)malloc(npixels*sizeof(struct carbonpool));  //49


	CNRold=(struct CNratio*)malloc(npixels*sizeof(struct CNratio));  //49
	CNRnew=(struct CNratio*)malloc(npixels*sizeof(struct CNratio));  //49


	soilw1_old= (float *)malloc(npixels*sizeof(float));//46  //第一层土壤水分含量 
	soilw2_old= (float *)malloc(npixels*sizeof(float));//47  //第二层土壤水分含量 
	soilw3_old= (float *)malloc(npixels*sizeof(float));//48  //第三层土壤水分含量 


	soilw1_n= (float *)malloc(npixels*sizeof(float));//46  //第一层土壤水分含量 
	soilw2_n= (float *)malloc(npixels*sizeof(float));//47  //第二层土壤水分含量 
	soilw3_n= (float *)malloc(npixels*sizeof(float));//48  //第三层土壤水分含量 



CI= (float *)malloc(npixels*sizeof(float));//49  //聚集度指数  
TI= (float *)malloc(npixels*sizeof(float));//50  //年平均温度



	for(i=0;i<npixels;i++) {

	lat[i]=0;    	lon[i]=0;      //经纬度中间变量, 短整型数
	
	lc[i]=0;     	awc[i]=0; 	lai[i]=0;   	laie[i]=0;
    
	Ndep0[i]=0; 	Ndep[i]=0;
	
	nppold[i]=0; 	npp[i]=0;    	nep[i]=0; 	hr[i]=0;
	//******lixuansong***** 
	output[i]=0; 	outputnep[i]=0; tnpp[i]=0; 	lambd[i]=0;thr[i]=0;lambd2[i]=0;

Cpoolold[i].Cstem=0;    Cpoolold[i].Cleaf=0; Cpoolold[i].Ccroot=0; Cpoolold[i].Cfroot=0;

Cpoolold[i].Ccd=0;      Cpoolold[i].Csmd=0;  Cpoolold[i].Cssd=0;   Cpoolold[i].Cfmd=0; 

 Cpoolold[i].Cfsd=0;    Cpoolold[i].Csm=0;   Cpoolold[i].Cm=0;    Cpoolold[i].Cs=0;     Cpoolold[i].Cp=0; 
	
CNRold[i].CNleaf=20;    CNRold[i].CNstem=20; CNRold[i].CNcroot=20;  CNRold[i].CNfroot=20;  CNRold[i].CNcd=20; CNRold[i].CNsmd=20; 

CNRold[i].CNssd=20;    CNRold[i].CNfmd=20;   CNRold[i].CNfsd=20;    CNRold[i].CNsm=20;     CNRold[i].CNm=20;  CNRold[i].CNs=20;	CNRold[i].CNp=20;  CNRold[i].Nav=20; 

CNRold[i].Nt=20; 



	soilw1_old[i]=0.1; 	soilw2_old[i]=0.1;   soilw3_old[i]=0.1; 
	}

	HY1=(struct hy_c1*)malloc(npixels*sizeof(struct hy_c1));  //49
	HY2=(struct hy_c2*)malloc(npixels*sizeof(struct hy_c2));  //50
	HY3=(struct hy_c3*)malloc(npixels*sizeof(struct hy_c3));   //51

	//lambdd=(float *)malloc(npixels*366*sizeof(float)); //环境因子影响土壤分解的系数
	//Tmean=(float *)malloc(npixels*366*sizeof(float));  //预热期间的日平均温度
	//NPPd=(float *)malloc(npixels*366*sizeof(float));  //预热期间的NPP
	
	//for(i=0;i<npixels;i++)lambdd[i]=(float *)malloc(366*sizeof(float)); 

	//*******	Daily Model Parameters  (Currently there are 24) ****
	x=(float*)malloc(SIZEX*sizeof(float));                          //52

	//*****	Memory for model parameters which change daily ****

	xx=(struct xvalue*)malloc(npixels*sizeof(struct xvalue));       //53

//sdat=(struct climatedata*)malloc(npixels*sizeof(struct climatedata));
    ClimateD=(short *) malloc(npixels5*sizeof(short));                 //54

	    
	    fseek(f_lat,1L*(long) line* (long) npixels*sizeof(float),SEEK_SET);  
	    fread(&lat[0],sizeof(float),npixels,f_lat);     
		
        fseek(f_lon,1L*(long) line* (long) npixels*sizeof(float),SEEK_SET); 
		fread(&lon[0],sizeof(float),npixels,f_lon);      
		

		for(pix=0; pix<npixels; pix++) {
		pixel = pix_offset + pix; 
		zeroxx1(pix,x,xx);              //交换中间变量 
		}	     


		printf("\n Firest SIMULATION IN LINE control=%d line=%d\n", control,line+1);	  

        printf("Land cover 1\n");

        f9=fopen("\\Global_Input\\Global_cover\\Global_0727_landcover.raw", "rb");
	    ptr=(line*npixels )*sizeof(char);   
        fseek(f9,ptr,SEEK_SET);
     	fread(&LC[0],1,npixels,f9);  
    	printf("Land cover 2\n");	
		for(pix=0; pix<npixels; pix++) lc[pix]=LC[pix];
    	fclose(f9);
        printf("Land cover3 \n");  
 
   
 
		f9=fopen("\\Global_Input\\Global_cover\\Global_CI_072727.img", "rb");
		ptr=(line*npixels )*sizeof(float);   
		fseek(f9,ptr,SEEK_SET);
		fread(&CI[0],4,npixels,f9);  
	    fclose(f9);
	

		f9=fopen("\\Global_Input\\Global_cover\\Global_Tmean_072727.dat", "rb");
		ptr=(line*npixels )*sizeof(float);   
		fseek(f9,ptr,SEEK_SET);
		fread(&TI[0],4,npixels,f9);  
		fclose(f9);



 //读入土壤属性数据  
		readsoildata(line,HY1,HY2,HY3,&rcode);              //
		if(rcode == ERROR) {
			printf("读入土壤属性数据出错");
			exit(0);
		}

  printf("Land cover4 \n");  

		for(pix=0; pix<npixels; pix++) { 
			// 初始化土壤含水量，假设土壤初始含水量等于田间持水量，该假设对干旱区而言，可能偏高

			soilw1_old[pix]= HY1[pix].FC;                    soilw2_old[pix]= HY2[pix].FC;        soilw3_old[pix]= HY3[pix].FC; 
			
			if(soilw1_old[pix]<0.025) soilw1_old[pix]=0.025; 	 if(soilw1_old[pix]>0.6) soilw1_old[pix]=0.6;
            if(soilw2_old[pix]<0.025) soilw2_old[pix]=0.025; 	 if(soilw2_old[pix]>0.6) soilw2_old[pix]=0.6;
			if(soilw3_old[pix]<0.025) soilw3_old[pix]=0.025; 	 if(soilw3_old[pix]>0.6) soilw3_old[pix]=0.6;
			
			//soilw2_old[pix]=__max(0.05,soilw2_old[pix]);     soilw2_old[pix]=__min(0.6,soilw2_old[pix]);
			//soilw3_old[pix]=__max(0.05,soilw3_old[pix]);     soilw3_old[pix]=__min(0.6,soilw3_old[pix]);

		}
		/////////////gpu==================================================================22222222222222222222222
		
		  printf("Land cover5 \n");  
		
		cudaMalloc((void**) &gpuHY1, npixels*sizeof(struct hy_c1));   //2
		cudaMalloc((void**) &gpuHY2, npixels*sizeof(struct hy_c2));   //3
		cudaMalloc((void**) &gpuHY3, npixels*sizeof(struct hy_c3));    //4 
		  printf("Land cover6 \n");  
		
		cudaMalloc((void**) &gpulat, npixels*sizeof(float));   //5
		cudaMalloc((void**) &gpulc, npixels*sizeof(int short)); //6
		
		cudaMalloc((void**) &gpusoilw1_old, npixels*sizeof(float)); //7
		cudaMalloc((void**) &gpusoilw2_old, npixels*sizeof(float)); //8
		cudaMalloc((void**) &gpusoilw3_old, npixels*sizeof(float)); //9


		cudaMalloc((void**) &gpusoilw1_n, npixels*sizeof(float)); //7
		cudaMalloc((void**) &gpusoilw2_n, npixels*sizeof(float)); //8
		cudaMalloc((void**) &gpusoilw3_n, npixels*sizeof(float)); //9


		cudaMalloc((void**) &gpuCI, npixels*sizeof(float)); //10
        cudaMalloc((void**) &gpuTI, npixels*sizeof(float)); //10

		cudaMemcpy(gpuHY1,HY1,npixels*sizeof(struct hy_c1),cudaMemcpyHostToDevice);//c2
		cudaMemcpy(gpuHY2,HY2,npixels*sizeof(struct hy_c2),cudaMemcpyHostToDevice);//c3
		cudaMemcpy(gpuHY3,HY3,npixels*sizeof(struct hy_c3),cudaMemcpyHostToDevice);//c4
		
		cudaMemcpy(gpulat,lat,npixels*sizeof(float),cudaMemcpyHostToDevice); //c5
		cudaMemcpy(gpulc,lc,npixels*sizeof(int short),cudaMemcpyHostToDevice);//c6
       	
		cudaMemcpy(gpusoilw1_old,soilw1_old,npixels*sizeof(float),cudaMemcpyHostToDevice);//c7
		cudaMemcpy(gpusoilw2_old,soilw2_old,npixels*sizeof(float),cudaMemcpyHostToDevice);//c8
		cudaMemcpy(gpusoilw3_old,soilw3_old,npixels*sizeof(float),cudaMemcpyHostToDevice);  ////c9

        cudaMemcpy(gpuCI,CI,npixels*sizeof(float),cudaMemcpyHostToDevice);//c10
        cudaMemcpy(gpuTI,TI,npixels*sizeof(float),cudaMemcpyHostToDevice);//c10
		
		
		cudaMalloc((void**) &gpuCpoolold,  npixels*sizeof(struct carbonpool));  //28

        cudaMalloc((void**) &gpuCpoolnew,  npixels*sizeof(struct carbonpool));  //28

		cudaMalloc((void**) &gpuCNRold,    npixels*sizeof(struct CNratio));  //28

		cudaMalloc((void**) &gpuCNRnew,    npixels*sizeof(struct CNratio));  //28

		
		
		
		
		
		
		
		
		for(pix=0; pix<npixels; pix++){
		xx[pix].x1=0;   	tnpp[pix]=0; 		lambd[pix]=0; thr[pix]=0; lambd2[pix]=0;
		   
	   }

		
		if(control==6){//***************************************************************************************************************************************
	
		cudaMalloc((void**) &gpulambd, npixels*sizeof(float));                                   //10
		cudaMemcpy(gpulambd,lambd,npixels*sizeof(float),cudaMemcpyHostToDevice);

		cudaMalloc((void**) &gpulambd2, npixels*sizeof(float));                                   //10
		cudaMemcpy(gpulambd2,lambd2,npixels*sizeof(float),cudaMemcpyHostToDevice);

       	//=======================================================================================第一轮年循环================================================================/
		for(year=startyear;year<=endyear1;year++){//年循环开始
		//===================================================================================================================================================================
			printf("\n year= %d...\n", year);	
			
				for(pix=0; pix<npixels; pix++) {
				pixel = pix_offset + pix; 
				npp[pix]=0.0; 	zeroxx(pix,x,xx);   output[pix]= 0;
				         
			} // end zfm赋初始值 

			//读入LAI年最大值用于估算生物量 
           tmpyear=startyear2;
        	readlaimax(line,xx,lc,&rcode,tmpyear);   

			if (rcode == ERROR) {
				printf ("读入最大LAI出错 "); exit(0);
			}
		
           		    
			co=c[0];   
			printf("%6.2f\n", co);

			//2013-10-25
			    cudaMalloc((void**) &gpusdat, npixels5*sizeof(short));    //11
				cudaMalloc((void**) &gpujday, sizeof(long));              //12 
				cudaMalloc((void**) &gpunpp, npixels*sizeof(float));      //14
			    cudaMalloc((void**) &gpuco, sizeof(float));                //16
			    cudaMalloc((void**) &gpulai, npixels*sizeof(float));       //18
     		    cudaMalloc((void**) &gpuxx, npixels*sizeof(struct xvalue));   //20
			    cudaMemcpy(gpuco, &co, sizeof(float), cudaMemcpyHostToDevice);
						
			for (jday=jday_start; jday<=jday_end; jday=jday++) {   //===================================================天时间循环开始=========================
	           if(jday%90==1) printf("jday=%d",jday);
				clock_t ttt1=clock();
                clock_t ttt2;//,ttt3,ttt4;

				//========================================读入气候数据============================================================================================		
			    //	getchar();

             readclim(year,jday,line,ClimateD,&rcode); //读入每天的气候数据
	
	          if (rcode == ERROR) {
		      printf("读入气候数据发生错误");	 exit(0);
              }


             ttt2=clock();
             if(jday==jday_start)
             printf("timh:%d\n",ttt2-ttt1);

			 
 //JUW2013  
     if(control1== 4 || control1==6){    // 1: all no change; 2: Climate change only; 3 :CO2 change only; 4 : LAI change only; 5:N_deposition,6: all change
	 
		 if(year<=2000){		 //2000年前读入AVHRR的LAI数据

	if( jday==1||jday==16||jday==32||jday==47||jday==60||jday==75||jday==91||jday==106||jday==121
		||jday==136||jday==152||jday==167||jday==182||jday==197||jday==213||jday==228||jday==244
		||jday==259||jday==274||jday==289||jday==305||jday==320||jday==335||jday==350)
	{	               
		tmpyear=year;              
		               
		//printf("tmpyear0=%d    \n",tmpyear);
		
		readlai(line,laie,jday,&rcode,tmpyear);     //读取有效叶面积指数
						if (rcode == ERROR){
							printf("读入LAI数据错误"); 	  
							exit(0);
 	
						}
						
						} 
   			}  //the end of year<2000


					else {  //2000年开始读入MODIS 的LAI数据
 
						if((jday-1)%8==0){ //判定是否读LAI数据，MODIS的LAI是8天一次 
							readlai(line,laie,jday,&rcode,tmpyear);     //读取有效叶面积指数
							if (rcode == ERROR){
								printf("读入LAI数据错误"); 	  
								exit(0);
							}

						}  // the end of if((jday-1)%8==0)
					
                      } //control1== 4    读入LAI数据结束
					
}  
    else{  //对于control==1,2,3的情景，
// 1: all no change; 2: Climate change only; 3 :CO2 change only; 4 : LAI change only; 5: all change
//========================================================================================================================================================================================//			
	
		tmpyear=year;   
		
		if( jday==1||jday==16||jday==32||jday==47||jday==60||jday==75||jday==91||jday==106||jday==121
			||jday==136||jday==152||jday==167||jday==182||jday==197||jday==213||jday==228||jday==244
			||jday==259||jday==274||jday==289||jday==305||jday==320||jday==335||jday==350)
	                 
	                    {
                       	readlai(line,laie,jday,&rcode,tmpyear);     //读取有效叶面积指数
						if (rcode == ERROR){
						printf("读入LAI数据错误"); 	  
						exit(0);
						} 
                    	}
//========================================================================================================================================================================================//
 	
}	 	
					////////gpu//lixuansong
					cudaMalloc((void**) &gpulaie, npixels*sizeof(float));      //21
					cudaMemcpy(gpulaie,laie,npixels*sizeof(float),cudaMemcpyHostToDevice);

					gpuLAI<<<BLOCK_NUM, THREAD_NUM,0>>>(gpulc,gpulaie,gpulai);
					cudaMemcpy(lai,gpulai,npixels*sizeof(float),cudaMemcpyDeviceToHost);
					cudaFree(gpulaie);  //right                     //21
   
                  //2013-10-25
			     
		
				cudaMemcpy(gpunpp,npp,npixels*sizeof(float),cudaMemcpyHostToDevice);
				cudaMemcpy(gpusdat,ClimateD,npixels*5*sizeof(short),cudaMemcpyHostToDevice);
				cudaMemcpy(gpujday,&jday,sizeof(long),cudaMemcpyHostToDevice);

				cudaMemcpy(gpuxx,xx,npixels*sizeof(struct xvalue),cudaMemcpyHostToDevice);
		      
               cudaMemcpy(gpusoilw1_old,soilw1_old,npixels*sizeof(float),cudaMemcpyHostToDevice);
               cudaMemcpy(gpusoilw2_old,soilw2_old,npixels*sizeof(float),cudaMemcpyHostToDevice);
               cudaMemcpy(gpusoilw3_old,soilw3_old,npixels*sizeof(float),cudaMemcpyHostToDevice);
			
gpuModel1<<<BLOCK_NUM, THREAD_NUM,0>>>(gpujday,gpupix_offset,gpulc,gpuCI,gpuTI,gpusdat,gpulai,gpunpp,gpuxx,gpulat,gpusoilw1_old,gpusoilw2_old,gpusoilw3_old,
									   gpuHY1,gpuHY2,gpuHY3,gpuco,gpulambd,gpulambd2);
		
				cudaMemcpy(npp,gpunpp,npixels*sizeof(float),cudaMemcpyDeviceToHost);
				cudaMemcpy(xx,gpuxx,npixels*sizeof(struct xvalue),cudaMemcpyDeviceToHost);
			           				 
				cudaMemcpy(soilw1_old,gpusoilw1_old,npixels*sizeof(float),cudaMemcpyDeviceToHost);
				cudaMemcpy(soilw2_old,gpusoilw2_old,npixels*sizeof(float),cudaMemcpyDeviceToHost);
				cudaMemcpy(soilw3_old,gpusoilw3_old,npixels*sizeof(float),cudaMemcpyDeviceToHost);
					 				
				j=j++;

                if(jday==jday_start)
                printf("time:%d\n",clock()-ttt2);
			
                 }//===========================//day循环结束===================================================================================

           //2013-10-25
			cudaFree(gpujday);  //13
			cudaFree(gpusdat);  //12
			cudaFree(gpunpp);   //14
						
			cudaFree(gpuxx);     //20
			cudaFree(gpulai);    //18
			cudaFree(gpuco);     //16
	
			for(pix=pix_offset;pix<npixels;pix++){
	    	tnpp[pix]=(float) tnpp[pix]+npp[pix] /(endyear1-startyear+1);
         	}
           



   }//=================================年循环结束===============================================================================
	 

				
		clock_t t0=clock();

	     
		cudaMalloc((void**) &gputnpp, npixels*sizeof(float));     //26
		cudaMemcpy(gputnpp,tnpp,npixels*sizeof(float),cudaMemcpyHostToDevice);
		
		
		//clock_t* time;
		
		//初始化土壤碳库
		 
		 gpuInit<<<BLOCK_NUM, THREAD_NUM,0>>>(gpupix_offset,gpuHY2,gpulc,gpuCI,gpuTI,gpulambd,gpulambd2,gputnpp,gpuCpoolold,gpuCNRold) ;
 
	 
	 	cudaMemcpy(Cpoolold,gpuCpoolold,npixels*sizeof(struct carbonpool),cudaMemcpyDeviceToHost);
		cudaMemcpy(CNRold,  gpuCNRold,  npixels*sizeof(struct CNratio),   cudaMemcpyDeviceToHost);
		 


        cudaFree(gpulambd);  //10
        cudaFree(gpulambd2);  //10

        cudaFree(gputnpp);    //26

	pix=1533;
		
		printf("Nt=%f,Nav=%f\n",CNRold[pix].Nt,CNRold[pix].Nav); 

		printf("\n");printf("\n");

	
		}  //the end of control==6====================================================================================================================================
  
else {  //对于其他情景，直接读入读入初始化数据
	  
	
	ptr=(lin*npixels )*sizeof(float); 
      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_stem%d_%d.img",startyear2-1,365); //1
	  f9=fopen(outfnpp_name,"rb");

      fseek(f9,ptr,SEEK_SET);
     
	  fread(&value[0],sizeof(float),npixels,f9);
      for(i=0;i<npixels;i++) Cpoolold[i].Cstem=value[i];
	  fclose(f9);  

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_leaf%d_%d.img",startyear2-1,365); //2
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
     
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++) Cpoolold[i].Cleaf=value[i];
	  fclose(f9);  

	

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_Ccroot%d_%d.img",startyear2-1,365); //3
	  f9=fopen(outfnpp_name,"rb");
     fseek(f9,ptr,SEEK_SET);

	 fread(&value[0],sizeof(float),npixels,f9);
	 for(i=0;i<npixels;i++) Cpoolold[i].Ccroot=value[i];
	 fclose(f9);  
  

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_froot%d_%d.img",startyear2-1,365); //4
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++) Cpoolold[i].Cfroot=value[i];
	  fclose(f9);  


      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_Cd%d_%d.img",  startyear2-1,365);  //5
	  f9=fopen(outfnpp_name,"rb");
     fseek(f9,ptr,SEEK_SET);
	 fread(&value[0],sizeof(float),npixels,f9);
	 for(i=0;i<npixels;i++) Cpoolold[i].Ccd=value[i];
	 fclose(f9);  
	 
	
      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_smd%d_%d.img", startyear2-1,365);  //6
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);

	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++) Cpoolold[i].Csmd=value[i];
	  fclose(f9);  

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_ssd%d_%d.img", startyear2-1,365);  //7
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++) Cpoolold[i].Cssd=value[i];
	  fclose(f9);  

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_fmd%d_%d.img",startyear2-1,365);  //8
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++) Cpoolold[i].Cfmd=value[i];
	  fclose(f9);  

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_fsd%d_%d.img",startyear2-1,365);  //9
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++) Cpoolold[i].Cfsd=value[i];
	  fclose(f9);  

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_sm%d_%d.img",startyear2-1,365);  //10
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++) Cpoolold[i].Csm=value[i];
	  fclose(f9);  

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_m%d_%d.img",startyear2-1,365);  //11
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++) Cpoolold[i].Cm=value[i];
	  fclose(f9);  

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_s%d_%d.img",startyear2-1,365);  //12
	  f9=fopen(outfnpp_name,"rb");
     fseek(f9,ptr,SEEK_SET);
	 fread(&value[0],sizeof(float),npixels,f9);
	 for(i=0;i<npixels;i++) Cpoolold[i].Cs=value[i];
	 fclose(f9);    

      sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_p%d_%d.img",startyear2-1,365);  //13
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++) Cpoolold[i].Cp=value[i];
	  fclose(f9);  

//==========================================================for C:N files========================================

      ptr=(lin*npixels )*sizeof(float); 
      
	  sprintf(outfnpp_name,"h:\\P7_output_all\\CN_stem%d_%d.img",startyear2-1,365); //1
	 
	
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
     

	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNstem=value[i];
	  fclose(f9);  

     sprintf(outfnpp_name,"h:\\P7_output_all\\CN_leaf%d_%d.img",startyear2-1,365); //2
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
     
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNleaf=value[i];
	  fclose(f9);  




      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_Ccroot%d_%d.img",startyear2-1,365); //3
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNcroot=value[i];
	  fclose(f9);  


      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_froot%d_%d.img",startyear2-1,365); //4
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNfroot=value[i];
	  fclose(f9);  


      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_Cd%d_%d.img",startyear2-1,365);  //5
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNcd=value[i];
	  fclose(f9);  


      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_smd%d_%d.img",startyear2-1,365);  //6
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNsmd=value[i];
	  fclose(f9);  


      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_ssd%d_%d.img",startyear2-1,365);  //7
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNssd=value[i];
	  fclose(f9);  


      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_fmd%d_%d.img",startyear2-1,365);  //8
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNfmd=value[i];
	  fclose(f9);  

      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_fsd%d_%d.img",startyear2-1,365);  //9
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNfsd=value[i];
	  fclose(f9);  


      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_sm%d_%d.img",startyear2-1,365);  //10
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNsm=value[i];
	  fclose(f9);  



      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_m%d_%d.img",startyear2-1,365);  //11
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNm=value[i];
	  fclose(f9);  



      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_s%d_%d.img",startyear2-1,365);  //12
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNs=value[i];
	  fclose(f9);  



      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_p%d_%d.img",startyear2-1,365);  //13
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].CNp=value[i];
	  fclose(f9);  



  
      sprintf(outfnpp_name,"h:\\P7_output_all\\Nav_%d_%d.img",startyear2-1,365);  //14
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
	  fread(&value[0],sizeof(float),npixels,f9);
	  for(i=0;i<npixels;i++)CNRold[i].Nav=value[i];
	  fclose(f9);  


 
      sprintf(outfnpp_name,"h:\\P7_output_all\\SW1_%d_%d.img",startyear2-1,365);  //15
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
      fread(&soilw1_old[0],sizeof(float),npixels,f9);
      fclose(f9);  

      
      sprintf(outfnpp_name,"h:\\P7_output_all\\SW2_%d_%d.img",startyear2-1,365);  //15
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
      fread(&soilw2_old[0],sizeof(float),npixels,f9);
      fclose(f9);  
  
      sprintf(outfnpp_name,"h:\\P7_output_all\\SW3_%d_%d.img",startyear2-1,365);  //15
	  f9=fopen(outfnpp_name,"rb");
      fseek(f9,ptr,SEEK_SET);
      fread(&soilw3_old[0],sizeof(float),npixels,f9);
      fclose(f9);  
     

}   //the end of control!=6 


    	//==================================================================================2018-03-07将下面几行关闭 
	//	for(pix=0; pix<npixels; pix++) {
		//	pixel = pix_offset + pix; 
			//zeroxx1(pix,x,xx);              //子程序
		//}	     






//===============================================================================================================================================
		printf("\n Second SIMULATION IN control=%d   LINE=%d...\n", control,line+1);	  

		// 	  Read landcover values (for each pixel in the line)  
		//readlc(line,lc,&rcode);                               
		//if (rcode == ERROR) {
		//	exit(0);
		//}

        f9=fopen("\\Global_Input\\Global_cover\\Global_0727_landcover.raw", "rb");
	    ptr=(line*npixels )*sizeof(char);   
        fseek(f9,ptr,SEEK_SET);
     	fread(&LC[0],1,npixels,f9);  for(pix=0; pix<npixels; pix++) lc[pix]=LC[pix];
        fclose(f9);
 

		f9=fopen("\\Global_Input\\Global_cover\\Global_CI_072727.img", "rb");
		ptr=(line*npixels )*sizeof(float);   
		fseek(f9,ptr,SEEK_SET);
		fread(&CI[0],4,npixels,f9);  
		fclose(f9);

		f9=fopen("\\Global_Input\\Global_cover\\Global_Tmean_072727.dat", "rb");
		ptr=(line*npixels )*sizeof(float);   
		fseek(f9,ptr,SEEK_SET);
		fread(&TI[0],4,npixels,f9);  
		fclose(f9);




       cudaMemcpy(gpulc,lc,npixels*sizeof(int short),cudaMemcpyHostToDevice);//c6

	   cudaMemcpy(gpuCI,CI,npixels*sizeof(float),cudaMemcpyHostToDevice);//c6
       cudaMemcpy(gpuTI,TI,npixels*sizeof(float),cudaMemcpyHostToDevice);//c6


    	// Read available water holding capacity for pixels in a line 
		readsoildata(line,HY1,HY2,HY3,&rcode); //            
		if(rcode == ERROR) {
		exit(0);
		}
		
		
		for(pix=0; pix<npixels; pix++){
			xx[pix].x1=0;   
				}

 /*	
f8=fopen("test.txt","wt");
fprintf(f8,"%s,%s,%s,%s,%s,","jday","Cstem","Ccroot","Cfroot", "Cleaf"); 
fprintf(f8,"%s,%s,%s,%s,%s,%s,%s,%s,%s,", "Ccd","Csmd","Cssd","Cfmd","Cfsd","Csm","Cm","Cs","Cp"); 
fprintf(f8,"%s,%s,%s,%s,","CNstem","CNcroot","CNfroot", "CNleaf"); 
fprintf(f8,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", "CNcd","CNsmd","CNssd","CNfmd","CNfsd","CNsm","CNm","CNs","CNp","Nav","NPP","Hr","NEP"); 
*/
 	

//==================================================开始正式模拟计算===============================================================================/

		if(control==6){//其他情景，直接读入数据
	
			for(pix=0; pix<npixels; pix++) { 
			// 初始化土壤含 水量，假设土壤初始含水量等于田间持水量，该假设对干旱区而言，可能偏高
			soilw1_old[pix]= HY1[pix].FC;                   soilw2_old[pix]= HY2[pix].FC;        soilw3_old[pix]= HY3[pix].FC; 
			//soilw1_old[pix]=__max(0.05,soilw1_old[pix]); 	soilw1_old[pix]=__min(0.6,soilw1_old[pix]);
			//soilw2_old[pix]=__max(0.05,soilw2_old[pix]);     soilw2_old[pix]=__min(0.6,soilw2_old[pix]);
			//soilw3_old[pix]=__max(0.05,soilw3_old[pix]);     soilw3_old[pix]=__min(0.6,soilw3_old[pix]);


			if(soilw1_old[pix]<0.025) soilw1_old[pix]=0.025; 	 if(soilw1_old[pix]>0.6) soilw1_old[pix]=0.6;
			if(soilw2_old[pix]<0.025) soilw2_old[pix]=0.025; 	 if(soilw2_old[pix]>0.6) soilw2_old[pix]=0.6;
			if(soilw3_old[pix]<0.025) soilw3_old[pix]=0.025; 	 if(soilw3_old[pix]>0.6) soilw3_old[pix]=0.6;

		}
	 }  //情景6结束

		

		cudaMemcpy(gpusoilw1_old,soilw1_old,npixels*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(gpusoilw2_old,soilw2_old,npixels*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(gpusoilw3_old,soilw3_old,npixels*sizeof(float),cudaMemcpyHostToDevice);
		//above================================== juw===========================================================================================================
			
	
		
		cudaMalloc((void**) &gpulambd, npixels*sizeof(float));    //55       
		cudaMemcpy(gpulambd,lambd,npixels*sizeof(float),cudaMemcpyHostToDevice);
	

    
		cudaMalloc((void**) &gpulambd2, npixels*sizeof(float));    //55       
		cudaMemcpy(gpulambd2,lambd2,npixels*sizeof(float),cudaMemcpyHostToDevice);


	
		srand(32730);
		for(year=1970;year<=endyear2;year++){//===================================================2013-10-06
        yearindex[year-startyear]=1970+rand()%10; 
       	}
	
       for(year=startyear;year<1980;year++){
        yearindex[year-startyear]=year;
        }


              cudaMalloc((void**) &gpunpp, npixels*sizeof(float));	  //15	
			  cudaMalloc((void**) &gpuhr, npixels*sizeof(float));       //23
			  cudaMalloc((void**) &gpunep, npixels*sizeof(float));       //24   
 		      cudaMalloc((void**) &gpuxx, npixels*sizeof(struct xvalue));    //21
		
		
		        cudaMalloc((void**) &gpuoutputnep, npixels3*sizeof( int short));  //25
                cudaMalloc((void**) &gpusdat, npixels5*sizeof(short));  //12
				cudaMalloc((void**) &gpujday, sizeof(long));            //13         
			


			//============================================================================================================================================
            ff2=fopen( "\\Global_Input\\Global_Ndeposition\\1891_I_Annual_dep_00727.dat","rb");
			ptr=sizeof( float)*(line*npixels + pix_offset);
			fseek(ff2,ptr,0);     
            fread(&Ndep0[0],4,npixels,ff2);
            fclose(ff2);
            if (rcode == ERROR) {
			exit(0);
			}


			for(pix=0; pix<npixels; pix++) {
				tnpp[pix]=0.0;  thr[pix]=0.0;
			}	
		//===================================================================================================================================================/

		//=========================================================================================第二轮年循环开始===========================================

		//===================================================================================================================================================/

		if(control==6) 	startyearII=startyear1;
		else 	        startyearII=startyear2;
			
		for(year=startyearII;year<=endyear2;year++){       //这里 startyear1=1871 
		printf("year=%d\n",year);
	
		for(pix=0; pix<npixels; pix++) {
		npp[pix]=0.0;  nppold[pix]=0.0;		nep[pix]=0.0;			hr[pix]=0.0; 		output[pix]= 0;
		zeroxx(pix,x,xx); 	
					}	                            // end zfm赋初始值 
        
		cudaMalloc((void**) &gpuNdep0, npixels*sizeof(float));    //56
		cudaMemcpy(gpuNdep0,Ndep0,npixels*sizeof(float),cudaMemcpyHostToDevice);


			//=======读入LAI maximum 数据，用于估算生物量=================================================================================================== 
			
             if(year<startyear2)  tmpyear=startyear2; 
			 else {
			if(control==4 || control==6 ) tmpyear=year;           //1: ALL no; 2: Climate change only; 3：CO2 change only;4: LAI 5 Ndeposion 6: all change; 
			else tmpyear=startyear2;   //JUW2013
			 }
			readlaimax(line,xx,lc,&rcode,tmpyear);    //读LAI年最大值
			if (rcode == ERROR) {
			exit(0);
			}
			
            //============================================================读入氮沉降数据======================================================================
			if(year<startyear2)  tmpyear=year; 
			 else {
			if(control==5 || control==6 ) tmpyear=year;           //1: ALL no; 2: Climate change only; 3：CO2 change only;4: LAI 5 Ndeposion 6: all change; 
			else tmpyear=startyear2;   //JUW2013
			 }

		    sprintf(filename, "\\Global_Input\\Global_Ndeposition\\%d%s",tmpyear,"_I_Annual_dep_00727.dat");
			ff2=fopen( filename,"rb");
			ptr=sizeof( float)*(line*npixels + pix_offset);
			fseek(ff2,ptr,0);     
            fread(&Ndep[0],4,npixels,ff2);
            fclose(ff2);
           //========================================================================================================================================================= 
		
			                  
            cudaMalloc((void**) &gpuNdep, npixels*sizeof(float));  //57
			cudaMemcpy(gpuNdep,Ndep,npixels*sizeof(float),cudaMemcpyHostToDevice);

		  if(year>=startyear2){                                        //startyear=1951 ;startyear0=1901
           if(control==3 || control==6) co=c[year-startyear0];           //1: ALL no; 2: Climate change only; 3：CO2 change only;4: LAI 5:N deposition 6: all change; 
		   else            	            co=c[startyear2-startyear0];      //co=c[startyear-startyear];   //2013-11-06
			}	
	    	else                        co=c[year-startyear0];  
	
//printf("year=%d co=%f %f %f\n",year,co,c[0],c[1]);


		    cudaMalloc((void**) &gpuco, sizeof(float));      //17
			cudaMemcpy(gpuco, &co, sizeof(float), cudaMemcpyHostToDevice);

                

//==================================================================================================================================================
 
			
for (jday=jday_start; jday<=jday_end; jday=jday++) {   //===============================时间循环开始=============================
                 if(jday%90==1) printf("jday=%d",jday);		

				 //start to read climate
	    		if (year>=startyear2){      //startyear2=1981
                if(control==2 || control==6  ) readclim(year,jday,line,ClimateD,&rcode);         //1: ALL no; 2: Climate change only; 3：CO2 change only;4: LAI 5 Ndepostion 6: all change; 
	                      else                 readclim(yearindex[year-startyear],jday,line,ClimateD,&rcode);			
		 
				}    //startyear2年以后
				else   {  //1980年以前
                
				  if(year<startyear)  year1=startyear+(year-startyear1)%10;	 //startyear:气候数据开始年			
				   else   year1=year;	

				   readclim(year1,jday,line,ClimateD,&rcode);                            
				}

						 
					
					if (rcode == ERROR) {
						printf("读入气候数据发生错误") ;	    
						exit(0);
					}
						
	        //the end of reading climate data			
 
                    
					if(control==4  || control==6	){//1: ALL no; 2: Climate change only; 3：CO2 change only;4: LAI 5 Ndepostion 6: all change; 
				if(year<=2000){
                  if(year<startyear2)      tmpyear=1980;                           //
				  else                     tmpyear=year;
			
				  if( jday==1||jday==16||jday==32||jday==47||jday==60||jday==75||jday==91||jday==106||jday==121
					  ||jday==136||jday==152||jday==167||jday==182||jday==197||jday==213||jday==228||jday==244
					  ||jday==259||jday==274||jday==289||jday==305||jday==320||jday==335||jday==350)

					readlai(line,laie,jday,&rcode,tmpyear);     //读取有效叶面积指数
					
									
					if (rcode == ERROR){
						printf("读入LAI数据错误"); 	  
						exit(0);
					} 
					}  //the end of year<2000

			else {
			if((jday-1)%8==0){ //判定是否读LAI数据，MODIS的LAI是8天一次
					
				tmpyear=year;
				readlai(line,laie,jday,&rcode,tmpyear);     //读取有效叶面积指数
						if (rcode == ERROR){
							printf("读入LAI数据错误"); 	  
							exit(0);

						}
					}  //the end of (jday-1)%8==0
					
			         }  // the end of elso : year>=2000
				 
					}   //the end of control==4 || control==5

 
					else{  //if(control==1 || control==2 || control==3 || control==5)

						tmpyear=1980;    
					
						if( jday==1||jday==16||jday==32||jday==47||jday==60||jday==75||jday==91||jday==106||jday==121
							||jday==136||jday==152||jday==167||jday==182||jday==197||jday==213||jday==228||jday==244
							||jday==259||jday==274||jday==289||jday==305||jday==320||jday==335||jday==350)
		  
		  
		  {
                     								
				  readlai(line,laie,jday,&rcode,tmpyear);     //读取有效叶面积指数
						if (rcode == ERROR){
							printf("读入LAI数据错误"); 	  
							exit(0);
						} 
                    	}  // the end of  if(control==1 || control==2 || control==3|| control==5 )
//========================================================================================================================================================================================//
            
					}
					
			cudaMalloc((void**) &gpulaie, npixels*sizeof(float));                              //22
			cudaMemcpy(gpulaie,laie,npixels*sizeof(float),cudaMemcpyHostToDevice);

			cudaMalloc((void**) &gpulai, npixels*sizeof(float));                                 //19
			gpuLAI<<<BLOCK_NUM, THREAD_NUM,0>>>(gpulc,gpulaie,gpulai);
			cudaMemcpy(lai,gpulai,npixels*sizeof(float),cudaMemcpyDeviceToHost);
			cudaFree(gpulaie);                                                                   //22
           
				
				cudaMemcpy(gpusdat,ClimateD,npixels5*sizeof(short),cudaMemcpyHostToDevice);
				cudaMemcpy(gpujday,&jday,sizeof(long),cudaMemcpyHostToDevice);
				cudaMemcpy(gpuxx,xx,npixels*sizeof(struct xvalue),cudaMemcpyHostToDevice);
 			
				
                cudaMemcpy(gpunpp,npp,npixels*sizeof(float),cudaMemcpyHostToDevice);
				cudaMemcpy(gpuhr,hr,npixels*sizeof(float),cudaMemcpyHostToDevice);
				cudaMemcpy(gpunep,nep,npixels*sizeof(float),cudaMemcpyHostToDevice);



				 if(jday==1){
	 
                cudaMemcpy(gpusoilw1_old,soilw1_old,npixels*sizeof(float),cudaMemcpyHostToDevice);
				cudaMemcpy(gpusoilw2_old,soilw2_old,npixels*sizeof(float),cudaMemcpyHostToDevice);
				cudaMemcpy(gpusoilw3_old,soilw3_old,npixels*sizeof(float),cudaMemcpyHostToDevice);		

				cudaMemcpy(gpuCpoolold,Cpoolold,npixels*sizeof(struct carbonpool),cudaMemcpyHostToDevice);
				
	            cudaMemcpy(gpuCNRold,CNRold,npixels*sizeof(struct CNratio),cudaMemcpyHostToDevice);

				  }

				 
				 
				 
				 
				 /*
	pix=1533;					
	printf("line=%d 1=%d 2=%d 3=%d 4=%d 5=%d \n",line,ClimateD[pix],ClimateD[pix+4950],ClimateD[pix+4950*2],ClimateD[pix+4950*3],ClimateD[pix+4950*4]);
    printf("line=%d 1=%f 2=%f 3=%f 4=%f 5=%f 6=%f\n",line,soilw1_old[pix],soilw2_old[pix],soilw1_old[pix],lai[pix],Ndep[pix],Ndep[pix]);   

 
	printf("year=%d,stem=%f,croot=%f,frrot=%f,leaf=%f\n",year,Cpoolold[pix].Cstem,Cpoolold[pix].Ccroot,Cpoolold[pix].Cfroot,Cpoolold[pix].Cleaf); 

	printf("cd=%f,smd=%f,ssd=%f,fmd=%f,fsd=%f,sm=%f,m=%f,s=%f,p=%f\n",Cpoolold[pix].Ccd,Cpoolold[pix].Csmd,Cpoolold[pix].Cssd,Cpoolold[pix].Cfmd,
		                                                              Cpoolold[pix].Cfsd,Cpoolold[pix].Csm,Cpoolold[pix].Cm,Cpoolold[pix].Cs,Cpoolold[pix].Cp); 
	
	
	printf("Nstem=%f,Ncroot=%f,Nfroot=%f,Nleaf=%f\n",CNRold[pix].CNstem,CNRold[pix].CNcroot,CNRold[pix].CNfroot, CNRold[pix].CNleaf); 
	
	
	
	printf("Ncd=%f,Nsmd=%f,Nssd=%f,Nfmd=%f,Nfsd=%f,Nsm=%f,Nm=%f,Ns=%f,Nt=%f,Nav=%f,npp=%f,hr=%f,nep=%f\n",CNRold[pix].CNcd,CNRold[pix].CNsmd,CNRold[pix].CNssd,CNRold[pix].CNfmd,
		                                           CNRold[pix].CNfsd,CNRold[pix].CNsm,CNRold[pix].CNm,CNRold[pix].CNs,CNRold[pix].Nt,CNRold[pix].Nav,npp[pix],hr[pix],nep[pix]); 

printf("\n");printf("\n");
*/

				gpuModel2<<<BLOCK_NUM, THREAD_NUM,0>>>(gpujday,gpupix_offset,gpulc,gpuCI,gpuTI,gpusdat,gpuNdep,gpuNdep0,gpulai,gpunpp,gpuxx,gpulat,
				gpusoilw1_old,gpusoilw2_old,gpusoilw3_old,gpuHY1,gpuHY2,gpuHY3,gpuco,gpulambd,gpulambd2,gpusoilw1_n,gpusoilw2_n,gpusoilw3_n,
				gpuCpoolold,gpuCNRold,gpuCpoolnew,gpuCNRnew,gpuhr,gpunep,gpuoutputnep);

//===================================================================================================================================================/
               // cudaMemcpy(lambd,gpulambd,npixels*sizeof(float),cudaMemcpyDeviceToHost);
			  	cudaMemcpy(xx,gpuxx,npixels*sizeof(struct xvalue),cudaMemcpyDeviceToHost);
     
				cudaMemcpy(hr,gpuhr,npixels*sizeof(float),cudaMemcpyDeviceToHost);
				cudaMemcpy(nep,gpunep,npixels*sizeof(float),cudaMemcpyDeviceToHost);
				cudaMemcpy(npp,gpunpp,npixels*sizeof(float),cudaMemcpyDeviceToHost);
				cudaMemcpy(outputnep,gpuoutputnep,npixels3*sizeof( int short),cudaMemcpyDeviceToHost);	


				cudaMemcpy(gpusoilw1_old,gpusoilw1_n,npixels*sizeof(float),cudaMemcpyDeviceToDevice);
				cudaMemcpy(gpusoilw2_old,gpusoilw2_n,npixels*sizeof(float),cudaMemcpyDeviceToDevice);
				cudaMemcpy(gpusoilw3_old,gpusoilw3_n,npixels*sizeof(float),cudaMemcpyDeviceToDevice);

				cudaMemcpy(gpuCpoolold,gpuCpoolnew,npixels*sizeof(struct carbonpool),cudaMemcpyDeviceToDevice);
				cudaMemcpy(gpuCNRold,  gpuCNRnew,npixels*sizeof(struct CNratio),cudaMemcpyDeviceToDevice);
			


//===============================================================================================================================================/
				  if(jday==365){

if(year<(startyear1+20)){
					  for(pix=0; pix<npixels; pix++) {

						  thr[pix]=hr[pix];
						  tnpp[pix]=npp[pix];

					  }	                            // end zfm赋初始值 
}
			    cudaMemcpy(soilw1_old,gpusoilw1_old,npixels*sizeof(float),cudaMemcpyDeviceToHost);
				cudaMemcpy(soilw2_old,gpusoilw2_old,npixels*sizeof(float),cudaMemcpyDeviceToHost);
		    	cudaMemcpy(soilw3_old,gpusoilw3_old,npixels*sizeof(float),cudaMemcpyDeviceToHost);	
					

				cudaMemcpy(Cpoolold,gpuCpoolold,npixels*sizeof(struct carbonpool),cudaMemcpyDeviceToHost);
				cudaMemcpy(CNRold,  gpuCNRold,npixels*sizeof(struct CNratio),cudaMemcpyDeviceToHost);

  }
			/*
    pix=903;
	fprintf(f8,"%d,%f,%f,%f,%f,",year,Cstem[pix],Ccroot[pix],Cfroot[pix], Cleaf[pix]); 
	fprintf(f8,"%f,%f,%f,%f,%f,%f,%f,%f,%f,", Ccd[pix],Csmd[pix],Cssd[pix],Cfmd[pix],Cfsd[pix],	Csm[pix],Cm[pix],Cs[pix],Cp[pix]); 
    fprintf(f8,"%f,%f,%f,%f,",CNstem[pix],CNcroot[pix],CNfroot[pix], CNleaf[pix]); 
	fprintf(f8,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", CNcd[pix],CNsmd[pix],CNssd[pix],CNfmd[pix],CNfsd[pix],CNsm[pix],CNm[pix],CNs[pix],Nt[pix],Nav[pix],npp[pix],hr[pix],nep[pix]); 
 	 */
  
/*
    pix=903;
	printf("year=%d,stem=%f,croot=%f,frrot=%f,leaf=%f\n",year,Cstem[pix],Ccroot[pix],Cfroot[pix], Cleaf[pix]); 
	printf("cd=%f,smd=%f,ssd=%f,fmd=%f,fsd=%f,sm=%f,m=%f,s=%f,p=%f\n", Ccd[pix],Csmd[pix],Cssd[pix],Cfmd[pix],Cfsd[pix],	Csm[pix],Cm[pix],Cs[pix],Cp[pix]); 
    printf("Nstem=%f,Ncroot=%f,Nfroot=%f,Nleaf=%f\n",CNstem[pix],CNcroot[pix],CNfroot[pix], CNleaf[pix]); 
	printf("Ncd=%f,Nsmd=%f,Nssd=%f,Nfmd=%f,Nfsd=%f,Nsm=%f,Nm=%f,Ns=%f,Nt=%f,Nav=%f,npp=%f,hr=%f,nep=%f\n", CNcd[pix],CNsmd[pix],CNssd[pix],CNfmd[pix],CNfsd[pix],CNsm[pix],CNm[pix],CNs[pix],Nt[pix],Nav[pix],npp[pix],hr[pix],nep[pix]); 
 
  		




				pix=1533;					
				printf("line=%d 1=%d 2=%d 3=%d 4=%d 5=%d \n",line,ClimateD[pix],ClimateD[pix+4950],ClimateD[pix+4950*2],ClimateD[pix+4950*3],ClimateD[pix+4950*4]);
				printf("line=%d 1=%f 2=%f 3=%f 4=%f 5=%f 6=%f\n",line,soilw1_old[pix],soilw2_old[pix],soilw1_old[pix],lai[pix],Ndep[pix],Ndep[pix]);   


				printf("year=%d,stem=%f,croot=%f,frrot=%f,leaf=%f\n",year,Cpoolold[pix].Cstem,Cpoolold[pix].Ccroot,Cpoolold[pix].Cfroot,Cpoolold[pix].Cleaf); 

				printf("cd=%f,smd=%f,ssd=%f,fmd=%f,fsd=%f,sm=%f,m=%f,s=%f,p=%f\n",Cpoolold[pix].Ccd,Cpoolold[pix].Csmd,Cpoolold[pix].Cssd,Cpoolold[pix].Cfmd,
					Cpoolold[pix].Cfsd,Cpoolold[pix].Csm,Cpoolold[pix].Cm,Cpoolold[pix].Cs,Cpoolold[pix].Cp); 


				printf("Nstem=%f,Ncroot=%f,Nfroot=%f,Nleaf=%f\n",CNRold[pix].CNstem,CNRold[pix].CNcroot,CNRold[pix].CNfroot, CNRold[pix].CNleaf); 



				printf("Ncd=%f,Nsmd=%f,Nssd=%f,Nfmd=%f,Nfsd=%f,Nsm=%f,Nm=%f,Ns=%f,Nt=%f,Nav=%f,npp=%f,hr=%f,nep=%f\n",CNRold[pix].CNcd,CNRold[pix].CNsmd,CNRold[pix].CNssd,CNRold[pix].CNfmd,
					CNRold[pix].CNfsd,CNRold[pix].CNsm,CNRold[pix].CNm,CNRold[pix].CNs,CNRold[pix].Nt,CNRold[pix].Nav,npp[pix],hr[pix],nep[pix]); 

				printf("\n");printf("\n");


*/




//===============================================================================================================================================/
	
				
 			if(((jday%30)==1 && year>=startyear) || (jday==365 && year>=startyear1) ) {  
	
//===============================================================================================================================================/
				 
			//============================输出NEP========================================       
					if(output_NPP=='y' || output_NPP=='Y') {       //按象元输出NeP         
		
					
					if(control==1)	sprintf(outfnpp_name,"h:\\P7_output_all_no\\gnep%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					if(control==2)	sprintf(outfnpp_name,"h:\\P7_output_climate\\gnep%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
                    if(control==3)	sprintf(outfnpp_name,"h:\\P7_output_CO2\\gnep%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j 
                    if(control==4)	sprintf(outfnpp_name,"h:\\P7_output_LAI\\gnep%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					if(control==5)  sprintf(outfnpp_name,"h:\\P7_output_Ndep\\gnep%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
                    if(control==6)  sprintf(outfnpp_name,"h:\\P7_output_all\\gnep%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
                    

					//	sprintf(outfnpp_name,".\\output\\gnep%d_%d.img",year,jday);    
						//sprintf(outfnpp_name,"cnep%d_%d.img",year,jday); 

					if(lin==0) {
						if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL)
						{
							printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
						}
						fwrite(&outputnep[0],sizeof( int short),npixels3,outfilenpp);

						//printf("1=%f   ", outputnep[0]);
						fclose(outfilenpp);
					
					}  //the end of lin==0
					
					else  {
						if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)
						{
							printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
						}
						fwrite(&outputnep[0],sizeof( int short),npixels3,outfilenpp);

						//printf("1=%f   ", outputnep[0]);
						fclose(outfilenpp);
					
					}
					
					}  //  //输出NEP结束  
//output Carbon pools
			}

 //=================================================开始保存C：N文件============================================


					//if(year>=startyear && control==6&&(jday==365)){
                  
							if( control==6&&(jday==365)){	
						
						sprintf(outfnpp_name,"h:\\P7_output_all\\CN_leaf%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
						if(lin==0) {
							if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
								printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
							}
						}		   
						else  			
						{
							if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
							printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
							}
							
							}
for(i=0;i<npixels;i++) value[i]=CNRold[i].CNleaf;
						fwrite(&value[0],sizeof(float),npixels,outfilenpp);
                       	fclose(outfilenpp);





						}


	   //      if(year>=startyear&& control==6 &&jday==365){    //2018-03-02


                     if(control==6 &&jday==365){    //2018-03-02
//============================================2==============================================================================
                      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_stem%d_%d.img",year,jday);     						
					  
					  
					  
					  if(lin==0) {
						  if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }
					  }		   
					  else  			
					  {
						  if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }

							}

					  for(i=0;i<npixels;i++) value[i]=CNRold[i].CNstem;
					  fwrite(&value[0],sizeof(float),npixels,outfilenpp);
					  fclose(outfilenpp);


//============================================3==============================================================================
                      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_Ccroot%d_%d.img",year,jday);    
					 
					  if(lin==0) {
						  if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }
					  }		   
					  else  			
					  {
						  if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }

							}


					  for(i=0;i<npixels;i++) value[i]=CNRold[i].CNcroot;
					  fwrite(&value[0],sizeof(float),npixels,outfilenpp);
					  fclose(outfilenpp);

				   
	
//============================================4==============================================================================
                      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_froot%d_%d.img",year,jday);    
					
					  if(lin==0) {
						  if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }
					  }		   
					  else  			
					  {
						  if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }

							}

					  for(i=0;i<npixels;i++) value[i]=CNRold[i].CNfroot;
					  fwrite(&value[0],sizeof(float),npixels,outfilenpp);
					  fclose(outfilenpp);

			
	
//============================================5==============================================================================
                     sprintf(outfnpp_name,"h:\\P7_output_all\\CN_Cd%d_%d.img",year,jday);       
					 
					 
					 if(lin==0) {
						 if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
							 printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							 getchar();
							 exit(0);
							}
					 }		   
					 else  			
						{
							if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
								printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
								getchar();
								exit(0);
							}

					 }


					 for(i=0;i<npixels;i++) value[i]=CNRold[i].CNcd;
					 fwrite(&value[0],sizeof(float),npixels,outfilenpp);
					 fclose(outfilenpp);
						
//============================================6=============================================================================
                  sprintf(outfnpp_name,"h:\\P7_output_all\\CN_smd%d_%d.img",year,jday);     
				 
				  if(lin==0) {
					  if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
						  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
						  getchar();
						  exit(0);
							}
				  }		   
				  else  			
						{
							if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
								printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
								getchar();
								exit(0);
							}

				  }

				  for(i=0;i<npixels;i++) value[i]=CNRold[i].CNsmd;
				  fwrite(&value[0],sizeof(float),npixels,outfilenpp);
				  fclose(outfilenpp);


//============================================7==============================================================================
                     sprintf(outfnpp_name,"h:\\P7_output_all\\CN_ssd%d_%d.img",year,jday);      
					 if(lin==0) {
						 if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
							 printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							 getchar();
							 exit(0);
							}
					 }		   
					 else  			
						{
							if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
								printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
								getchar();
								exit(0);
							}

					 }
					 for(i=0;i<npixels;i++) value[i]=CNRold[i].CNssd;
					 fwrite(&value[0],sizeof(float),npixels,outfilenpp);
					 fclose(outfilenpp);


						
//============================================8==============================================================================
                   sprintf(outfnpp_name,"h:\\P7_output_all\\CN_fmd%d_%d.img",year,jday);      
				   if(lin==0) {
					   if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
						   printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
						   getchar();
						   exit(0);
					   }
				   }		   
				   else  			
				   {
					   if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
						   printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
						   getchar();
						   exit(0);
					   }

							}

				   for(i=0;i<npixels;i++) value[i]=CNRold[i].CNfmd;
				   fwrite(&value[0],sizeof(float),npixels,outfilenpp);
				   fclose(outfilenpp);


					
//============================================9==============================================================================
                sprintf(outfnpp_name,"h:\\P7_output_all\\CN_fsd%d_%d.img",year,jday);      
				if(lin==0) {
					if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
						printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
						getchar();
						exit(0);
					}
				}		   
				else  			
				{
					if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
						printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
						getchar();
						exit(0);
					}

				}
				for(i=0;i<npixels;i++) value[i]=CNRold[i].CNfsd;
				fwrite(&value[0],sizeof(float),npixels,outfilenpp);
				fclose(outfilenpp);


				//============================================10=============================================================================
                sprintf(outfnpp_name,"h:\\P7_output_all\\CN_sm%d_%d.img",year,jday);      
				if(lin==0) {
					if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
						printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
						getchar();
						exit(0);
					}
				}		   
				else  			
				{
					if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
						printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
						getchar();
						exit(0);
					}

				}

				for(i=0;i<npixels;i++) value[i]=CNRold[i].CNsm;
				fwrite(&value[0],sizeof(float),npixels,outfilenpp);
				fclose(outfilenpp);

//============================================11============================================================================
            sprintf(outfnpp_name,"h:\\P7_output_all\\CN_m%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
			if(lin==0) {
				if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
					printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
					getchar();
					exit(0);
				}
			}		   
			else  			
			{
				if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
					printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
					getchar();
					exit(0);
				}

			}
			for(i=0;i<npixels;i++) value[i]=CNRold[i].CNm;
			fwrite(&value[0],sizeof(float),npixels,outfilenpp);
			fclose(outfilenpp);
			
		
//============================================12==============================================================================
                      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_s%d_%d.img",year,jday);  					
					  if(lin==0) {
						  if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }
					  }		   
					  else  			
					  {
						  if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }

							}
					  for(i=0;i<npixels;i++) value[i]=CNRold[i].CNs;
					  fwrite(&value[0],sizeof(float),npixels,outfilenpp);
					  fclose(outfilenpp);

//============================================13==============================================================================
                      sprintf(outfnpp_name,"h:\\P7_output_all\\CN_p%d_%d.img",year,jday);     
					  if(lin==0) {
						  if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }
					  }		   
					  else  			
					  {
						  if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }

							}
					  for(i=0;i<npixels;i++) value[i]=CNRold[i].CNp;
					  fwrite(&value[0],sizeof(float),npixels,outfilenpp);
					  fclose(outfilenpp);
					
//============================================14==============================================================================
                      sprintf(outfnpp_name,"h:\\P7_output_all\\Nav_%d_%d.img",year,jday);     
					  if(lin==0) {
						  if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }
					  }		   
					  else  			
					  {
						  if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
							  printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							  getchar();
							  exit(0);
						  }

							}

					  for(i=0;i<npixels;i++) value[i]=CNRold[i].Nav;
					  fwrite(&value[0],sizeof(float),npixels,outfilenpp);
					  fclose(outfilenpp);

					
	//============================================15==============================================================================JUW
                      sprintf(outfnpp_name,"h:\\P7_output_all\\SW1_%d_%d.img",year,jday);     
                       if(lin==0) {
						if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL)
						{
							printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
						}
                    fwrite(&soilw1_old[0],sizeof(float),npixels,outfilenpp);		
					
						
					  //  fwrite(&Nt[0],sizeof(float),npixels,outfilenpp);		
						
						fclose(outfilenpp);
					
					}  //the end of lin==0
					
					else  {
						if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)
						{
							printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
						}
                      fwrite(&soilw1_old[0],sizeof(float),npixels,outfilenpp);		
                      fclose(outfilenpp);
					}	

//============================================16==============================================================================
                      sprintf(outfnpp_name,"h:\\P7_output_all\\SW2_%d_%d.img",year,jday);     
                       if(lin==0) {
						if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL)
						{
							printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
						}
                    fwrite(&soilw2_old[0],sizeof(float),npixels,outfilenpp);		
					fclose(outfilenpp);
					
					}  //the end of lin==0
					
					else  {
						if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)
						{
							printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
						}
                      fwrite(&soilw2_old[0],sizeof(float),npixels,outfilenpp);		
                      fclose(outfilenpp);
					}	


//============================================17==============================================================================
                      sprintf(outfnpp_name,"h:\\P7_output_all\\SW3_%d_%d.img",year,jday);     
                       if(lin==0) {
						if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL)
						{
							printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
						}
                    fwrite(&soilw3_old[0],sizeof(float),npixels,outfilenpp);		
					fclose(outfilenpp);
					
					}  //the end of lin==0
					
					else  {
						if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)
						{
							printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
							getchar();
							exit(0);
						}
                      fwrite(&soilw3_old[0],sizeof(float),npixels,outfilenpp);		
                      fclose(outfilenpp);
					}	
                     }  //the end  if(control==6 && year==(startyear2-1))

//=======================================保存C：N文件结束=============================================

					 
//===========================================开始保存C库文件======================================================					 
					 if(control==6 &&jday==365){
//================================================================1================================================
                   sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_leaf%d_%d.img",year,jday);      
					
				   if(lin==0) {
					   if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
						   printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
						   getchar();
						   exit(0);
					   }
				   }		   
				   else  			
				   {
					   if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
						   printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
						   getchar();
						   exit(0);
					   }

							}
				   for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Cleaf;
				   fwrite(&value[0],sizeof(float),npixels,outfilenpp);
				   fclose(outfilenpp);

					
//============================================2==============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_stem%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}	 
                 
for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Cstem;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);

//============================================3==============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_Ccroot%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}

for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Ccroot;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);

//============================================4==============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_froot%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}
		
for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Cfroot;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);


//============================================5==============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_Cd%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}

for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Ccd;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);


//============================================6=============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_smd%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}
for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Csmd;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);					


//============================================7==============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_ssd%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}
	

for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Cssd;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);

//============================================8==============================================================================
sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_fmd%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j


if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}						
     
for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Cfmd;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);


//============================================9==============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_fsd%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}

for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Cfsd;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);


//============================================10=============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_sm%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}

for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Csm;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);

//============================================11============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_m%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}

for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Cm;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);

//============================================12==============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_s%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}

for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Cs;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);


//============================================13==============================================================================

sprintf(outfnpp_name,"h:\\P7_output_all\\Carbon_p%d_%d.img",year,jday);      ///////sprintf 连接"cnpp90_%d.img"+j
					
if(lin==0) {
	if ((outfilenpp=fopen(outfnpp_name, "wb"))== NULL){	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}
}		   
else  			
{
	if ((outfilenpp=fopen(outfnpp_name, "ab"))== NULL)	{	
		printf("\n Unable to open file <%s>,  exitting program ...\n\n", outf_name);
		getchar();
		exit(0);
	}

}

for(i=0;i<npixels;i++) value[i]=Cpoolold[i].Cp;
fwrite(&value[0],sizeof(float),npixels,outfilenpp);
fclose(outfilenpp);


//=====================================================================================================================/


}  //  //输出carbon pools结束  


					j=j++;
 


				
				
				cudaFree(gpulai);                              //19



	}    //=====================================================day循环结束=======================
			
     
			cudaFree(gpuco);  //17
            cudaFree(gpuNdep);  //57
            cudaFree(gpuNdep0);  //56


if(year<(startyear1+20)){

for(pix=0;pix<npixels;pix++){
if(lc[pix]!=0&& lc[pix]!=21&& lc[pix]!=19 ) {

	
	coef_c(pix,lc[pix],coef,HY2);
ratio=tnpp[pix]/thr[pix];
/*
if(pix==3666){
printf("tnpp=%f  hr=%f %f\n",tnpp[pix],thr[pix],ratio);

printf("%f  %f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f %f\n",

Cpoolold[pix].Cstem,Cpoolold[pix].Ccroot,Cpoolold[pix].Cfroot,Cpoolold[pix].Cleaf,Cpoolold[pix].Ccd,Cpoolold[pix].Csmd,Cpoolold[pix].Cssd,

Cpoolold[pix].Cfmd,Cpoolold[pix].Cfsd,Cpoolold[pix].Csm,Cpoolold[pix].Cm,Cpoolold[pix].Cs,Cpoolold[pix].Cp);

}
*/

/*
Cpoolold[pix].Cstem=Cpoolold[pix].Cstem*ratio;
Cpoolold[pix].Ccroot=Cpoolold[pix].Ccroot*ratio;

Cpoolold[pix].Cfroot=Cpoolold[pix].Cfroot*ratio;
Cpoolold[pix].Cleaf=Cpoolold[pix].Cleaf*ratio;
*/

Cpoolold[pix].Ccd=Cpoolold[pix].Ccd*ratio;

Cpoolold[pix].Csmd=Cpoolold[pix].Csmd*ratio;

Cpoolold[pix].Cssd=Cpoolold[pix].Cssd*ratio;

Cpoolold[pix].Cfmd=Cpoolold[pix].Cfmd*ratio;

Cpoolold[pix].Cfsd=Cpoolold[pix].Cfsd*ratio;

Cpoolold[pix].Csm=Cpoolold[pix].Csm*ratio;

Cpoolold[pix].Cm=Cpoolold[pix].Cm*ratio;

Cpoolold[pix].Cs=Cpoolold[pix].Cs*ratio;

Cpoolold[pix].Cp=Cpoolold[pix].Cp*ratio;

/*
if(pix==3666){
	printf("tnpp=%f  hr=%f %f\n",tnpp[pix],thr[pix],ratio);

	printf("%f  %f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f %f\n",

		Cpoolold[pix].Cstem,Cpoolold[pix].Ccroot,Cpoolold[pix].Cfroot,Cpoolold[pix].Cleaf,Cpoolold[pix].Ccd,Cpoolold[pix].Csmd,Cpoolold[pix].Cssd,

		Cpoolold[pix].Cfmd,Cpoolold[pix].Cfsd,Cpoolold[pix].Csm,Cpoolold[pix].Cm,Cpoolold[pix].Cs,Cpoolold[pix].Cp);


pix=pix;


}

*/


}
}// pix loop
}// year<(startyear1+20


if(year<(startyear1+10)){

	for(pix=0;pix<npixels;pix++){
	if(lc[pix]!=0&& lc[pix]!=21&& lc[pix]!=19 ) {
	coef_c(pix,lc[pix],coef,HY2);

   CNRold[pix].CNcd  =coef[35];     //Cd  库碳氮比
	CNRold[pix].CNssd =coef[36];     //ssd 库碳氮比
	CNRold[pix].CNsmd =coef[37];     //smd 库碳氮比 
	CNRold[pix].CNfsd =coef[38];     //fsd 库碳氮比
	CNRold[pix].CNfmd  =coef[39] ;    //fmd 库碳氮比
	CNRold[pix].CNsm  =coef[41];     //sm  库碳氮比
    CNRold[pix].CNm    =coef[42];     //sm  库碳氮比
   CNRold[pix]. CNs    =coef[40];     //slow库碳氮比
    CNRold[pix].CNp    =coef[43];     
    CNRold[pix].CNstem =coef[47];    //Setm C:NinitialValues[25] ;
	CNRold[pix].CNcroot=coef[47];   //Setm C:NinitialValues[25] ;
	CNRold[pix].CNleaf  =coef[46] ;   //叶C:N
	CNRold[pix].CNfroot =coef[46];   //叶C:N 
}
}
} // year<(startyear1+10

//JUW1


}//====================================年循环结束======================================================

          
           

            cudaFree(gpunpp);   //15
			cudaFree(gpuhr);    //23
			cudaFree(gpunep);   //24
			cudaFree(gpuxx);         //21




            cudaFree(gpuoutputnep);  //25
            cudaFree(gpujday);  //13
			cudaFree(gpusdat);  //12
			


			cudaFree(gpuCpoolold);
cudaFree(gpuCpoolnew);
				

cudaFree(gpuCNRold);
cudaFree(gpuCNRnew);
				cudaFree(gpusoilw1_n);//7
				cudaFree(gpusoilw2_n); //8
				cudaFree(gpusoilw3_n);  //9

  	           cudaFree(gpuHY1);   //2
		       cudaFree(gpuHY2);   //3
		       cudaFree(gpuHY3);   //4
		       cudaFree(gpulat);  //5
		       cudaFree(gpulc);   //6
		       cudaFree(gpusoilw1_old);//7
		       cudaFree(gpusoilw2_old); //8
    	       cudaFree(gpusoilw3_old);  //9
	   
 cudaFree(gpuCI);   //10
cudaFree(gpuTI);
		        cudaFree(gpulambd);   //55
  cudaFree(gpulambd2);   //55



	printf("time used: %d\n", clock()-ttt);
	free(LC);       printf ("1");  free(CI);   free(TI);  
	free(lc);       printf ("2");  
	free(awc);		printf ("3");
	free(lai);	    printf ("4");
	free(lon);	    printf ("5");
	free(lat);	    printf ("6");
	free(npp);	    printf ("7");
	
	free(output);	printf ("8");
	free(outputnep); printf ("9"); 
	free(ClimateD);  //10
 
     
     free(laie);     //12
	 free(nppold);   //13
	 free(nep);      //14
	 free(hr);       //15
 	 free(tnpp);     //16
    free(thr);     //16

	 free(lambd);    //17
 free(lambd2);    //17
	
	free(soilw1_old);//46  
	free(soilw2_old);   //47
	free(soilw3_old);   //48


	free(soilw1_n);//46  
	free(soilw2_n);   //47
	free(soilw3_n);   //48


free(Cpoolold);
free(Cpoolnew);

free(CNRold);
free(CNRnew);
free(value);


	free(HY1);  //49
	free(HY2); //50
	free(HY3); //51


 	free(x) ; //52
 	free(xx);  //53

    free(Ndep0) ; //54
 	free(Ndep);  //55

    //fclose(f8);
//exit(0);

 }  // ==================================================================行循环结束 zfm4.22==============================================================
cudaFree(gpupix_offset);  //1

} // the end of control loop

return 1;

}	/* end of main */





