/*************************************************************************
  program: 	beps.h
  Description:  Header file for defining constants and global variables
  ----------    for BEPS program

***************************************************************************
  CCRS (EMS/Applications Division)
  Wrintten by: 	J. Liu
  Modified by:  X. F. Feng
  Last update:  July 2003
*****************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<malloc.h>
#include <cuda_runtime.h>
#include <cutil_inline.h>
//#include <shrQATest.h>
#include <cuda.h>
/* Define Constants */

#define NOERROR		0
#define ERROR	    1
#define PI	    	 3.1415926

#define linoffset    0
#define pixoffset    0


#define nlines       5         //计算区域的行数 
#define npixels      2069100    //计算区域的列数
#define npixels3     6207300
#define npixels5     10345500


/*
#define nlines       5         //计算区域的行数 
#define npixels      2069100    //计算区域的列数
#define npixels3     6207300
#define npixels5     10345500




#define nlines       3         //计算区域的行数 
#define npixels      3450150   //计算区域的列数
#define npixels3     10350450
#define npixels5     17250750


#define nlines       5         //计算区域的行数 
#define npixels      2069100    //计算区域的列数
#define npixels3     6207300
#define npixels5     10345500


#define nlines       2091                      //计算区域的行数 
#define npixels      4950                      //计算区域的列数
#define npixels3     14850
#define npixels5     24750

#define nlines       418                    //计算区域的行数 
#define npixels      24750                     //计算区域的列数
#define npixels3     74250
#define npixels5     123750

#define nlines       209                      //计算区域的行数 
#define npixels      49500                      //计算区域的列数
#define npixels3     148500
#define npixels5     247500



#define nlines       10                   //计算区域的行数 
#define npixels      1034550//1291950       //计算区域的列数
#define npixels3     3103650//3875850
#define npixels5     5172750//6459750


*/


#define endday    15
#define factor0   10000

#define layer        3
#define jday_start   1    // Julian day start 

#define  jday_end  365  //lixuansong

#define startyear0  1871    //CO2数据开始年

#define startyear1  1871    //第二次正式计算的开始年  

#define startyear   1901    //第一次预热计算的起始年

#define endyear1    1905    //第一次预热计算的结束年

#define startyear2   1981   //LAI数据的开始年
#define endyear2     2016   //计算结束年 




#define depth1       0.15
#define depth2       0.45    //原为0.25
//#define depth3     改为在gpudoflux.cu子程序中随植被类型变化   

#define SIZEB		100
#define SIZEG 		100
#define SIZEX 		100
#define SIZEZ 		100

#define RTIMES 		1
#define RPERIOD 	3600	// in sec 


//bool InitCUDA();


void solution11(float xx[]);

void v_moisture(int long pix,  float B[],float K[],	float Psi[],struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[],float *theta_o,float *theta_n);

void inter(float  b[],float g[],float x[],float z[]);

void melt (float b[],float g[],float x[],float z[]);

void rs(int long pix,float b[],float g[],float x[],float z[],float rr[],struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[]);

float penmon(float tair,float vpd,float wrad,float acond,float scond);

void readxx(int long pix,float x[],struct xvalue xx[]);

void writexx(int long pix, float x[],struct xvalue xx[]);

void doflux(int long pix,float b[],float g[],float x[],float z[],float co,float soilw1,float soilw2,float soilw3,struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[],
	        float *soilwn1,float*soilwn2,float*soilwn3);


void gpulw1(float FC1,float FC2, float FC3, float sw1, float  sw2, float sw3, float *lw1,float *lw2, float *lw3);



void readb_init(float b[]); 

void zeroxx1(int long pix,float x[],struct xvalue xx[]);

void zeroxx(int long pix,float x[],struct xvalue xx[]);

void readlc(int long line,int short *lc,int short *rcode);

void readsoildata(long line,struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[],int short *rcode);

void readlaimax(int long line,struct xvalue xx[],int short *lc,int short *rcode,long tmpyear);

void readclim(long year,long jday,long line,short *climatedata,int short *rcode);

void readlai(int long line,float *lai,long jday,int short *rcode,long tmpyear);

void readb(float *b, short int lc_p);

void coef_c(long pix,short lc,float coef[], struct hy_c2 HY2[]);

void model(long jday,int long pix,float lat_p,float lai_p,short int lc_p,float CI,float TI,float sw_o1,float sw_o2,float sw_o3,struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[],
           float x[],float b[],short  climatedata[],struct  xvalue xx[],float co);

void readclim10(long year,long jday,long line,struct climatedata ****sdat,int short *rcode);

void preInitials(float  NPP,float lambd,float lambddd[],float coef1[],float initialValues[]);

void soilresp( float coef[],float  lambda,float Cstem1,float Cleaf1,float Ccroot1,float Cfroot1,float Ccd1,float Csmd1,float Cssd1,float Cfmd1 ,float Cfsd1,
              float Csm1 ,float Cm1 ,float Cs1 ,float Cp1,float  initialValues[]);

int rad_ssl(float sg,float cos_theta,float lai_p,float omega,float *ssun,float *sshade);

int net_shortwave(float sg,float cos_theta,float lai_p,float omega,float *ssun,float *sshade,float lai_u,float *n_ssun, float *n_sshade,float *n_sunder,float *n_sground);

int net_longwave (float lai_o, float lai_u,float omega,float es,float ta,float *m_lnet_o,float *m_lnet_u, float *lnet_g);		

void carbon(float b[],float g[],float x[],float z[],float co);

int farq_psn(float b[],float z[],float g[],float co);

//gpu
__device__  void gpumodel(float CNl,long jday,int long pix,long index,float lat_p,float lai_p, int short lc_p,float CI,float TI,float sw_o1,float sw_o2,float sw_o3,struct hy_c1 HY1[],
                           struct hy_c2 HY2[],struct hy_c3 HY3[], float x[], float  b[],short climatedata[],struct  xvalue xx[],
	                       float co,float*soilw1_new,float* soilw2_new,float* soilw3_new,float *tmean,float *t_d,float *t_n,float *d_l);

__device__ void gpusolution11(float xx[],float x1,float x2[],float x3[],float Y[]);

__device__ void gpucarbon(float CNl,long index,float b[],float g[],float x[],float z[],float co,float TI);

__device__ int gpufarq_psn(float CNl,float b[],float z[],float g[],float co);

__device__ void gpucoef_c(long pix, int short lc,float coef[],struct hy_c2 HY2[]);

__device__ void gpudoflux(float CNl,int long pix,long index,float b[],float g[],float x[],float z[],float co,float TI,float soilw1,float soilw2,float soilw3,struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[],float *soilwn1,float*soilwn2,float*soilwn3);

__device__  void gpulw1(float FC1,float FC2, float FC3, float sw1, float  sw2, float sw3, float *lw1,float *lw2, float *lw3);

__device__  void gpusetx( int short lc_p, float awc_p, float x[], float b[]);

__device__ void gpuzcomp(long jday,int long pix,float lat_p, int short lc_p,float CI,float TI,float awc_p,float lai_p,float b[],float   x[],float z[],short climatedata [],float*tmean);

__device__ void gpuinter(float  b[],float g[],float x[],float z[]);

__device__  void gpuzeroxx1(int long pix,float x[],struct xvalue xx[]);

__device__  void gpuzeroxx(int long pix,float x[],struct xvalue xx[]);

__device__  void gpureadxx( int long pix,  float x[],  struct xvalue xx[]);

__device__  void gpuwritexx(  int long pix,  float x[],  struct xvalue xx[]);

__device__ void gpumelt (float b[],float g[],float x[],float z[]);

__device__ float gpupenmon(float tair,float vpd,float wrad,float acond,float scond);

__device__ void   gpupreInitials(float NPP,float lambd,float lambddd[], float Tmean_m[], float NPP_d[],float coef1[], float initialValues[]);

__device__  void gpureadb(float *b, int short lc_p);

__device__ void gpurs(int long pix,float b[],float g[],float x[],float z[],float rr[],float *rrr1,float *rrr2,float *rrr3,float TI,struct hy_c1 HY1[],struct hy_c2 HY2[], struct hy_c3 HY3[]);

__device__ void gpuv_moisture(int long pix, float depth3, float B[],float K[],float Psi[],struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[],
                              float *theta_o,float *theta_n1,float *theta_n2,float *theta_n3);

__device__ int gpunet_shortwave(float sg, float cos_theta,float lai_p, int short lc_p,float omega,float *ssun,float *sshade,float lai_u,float *n_ssun,float *n_sshade,float *n_sunder,		 
	                            float *n_sground);
	
__device__ int gpunet_longwave(float lai_o,float lai_u,float omega,float es,float ta,float *m_lnet_o,float *m_lnet_u,float *lnet_g);

__device__ int gpurad_ssl(float sg,float cos_theta, float lai_p, int short lc_p, float omega,float *ssun,float *sshade);



__device__ void  gpusoilresp(int long pix,int short lc,float coef[],long jday,float lat,float dTmean,float dNPP,float Ndep, float Ndep0,
                             struct carbonpool Cpoolold[], struct CNratio CNRold[], struct carbonpool Cpoolnew[], struct CNratio CNRnew[],
                             float lambda, float lambd2); 


__device__ void gpureadb_init(float b[]);


// Declare structures 
/*
struct climatedata
	{
    short precip;
	short rehum;
	short srad;
	short tmax;
	short tmin;
  };
*/

short ClimateD;
struct hy_c1
	{  
	float PR;
	float FC;
	float WP;
	float SP;
	float B;
	float K; 
	float clay;
	float silt;
	};

struct hy_c2
	{  float PR; float FC;float WP;	float SP;float B;float K;  float clay; float silt; };

struct hy_c3
	{float PR; float FC;float WP;float SP;float B;	float K; float clay; float silt;   };


struct carbonpool
{
float Cstem;
float Cleaf;
float Ccroot;
float Cfroot;
float Ccd;
float Csmd;
float Cssd;
float Cfmd;
float Cfsd;
float Csm;
float Cm;
float Cs;
float Cp;
float Hr;
};


struct CNratio
{

float CNstem;
float CNleaf;
float CNcroot;
float CNfroot;
float CNcd;
float CNsmd;
float CNssd;
float CNfmd;
float CNfsd;
float CNsm;
float CNm;
float CNs;
float CNp;
float Nav;
float Nt;
};







struct xvalue
	{
	int long pix;
	float x1;      /* snow pack */
	float x2;		/* soil water */
	float x3;		/* out flow */
	float x4;		/* transpiration */
	float x5;		/* evaporation */
	float x6;		/* photosynthesis */
	float x7;		/* Evaporation from canopy */
	float x8;		/* leaf C */
	float x9;		/* stem C */
	float x10;		/* fine root C */
	float x11;		/* GPP */	
	float x12;		/* ResSoil */
	float x13;		/* NEP */
	float x14;		/* sublimation from canopy*/
	float x15;		/* Evaporation from soil  */
	float x16;		/* sublimation from soil  */
	float x17;     /* Number of days with rain */
	float x18;		/* understory Trans */
	float x19;     /*  canopy trans*/
	float x21;     /* total Res */
    float x22;     /* M Res leaf */
	float x23;     /* M Res stem */
	float x24;     /* M Res stem */
	float x25;     /* M Res stem */
	
	float x31;      /* leaf*/
	float x32;      /* stem*/
	float x33;      /* fine root*/
	float x34;      /* coarse root*/
	
	float x41;      /* leaf*/
	float x42;      /* stem*/
	float x43;      /* fine root*/
	float x44;      /* coarse root*/
	
	
};




/*	Declare global variables */

char awcf[255];         /* Available water holding capacity file */
char landcoverf[255];   /* Land cover file */
char leafareaf[255];    /* Leaf area index file */
char snowpackf[255];    /* Snow pack data file, XF, 24 June, 2003*/ 
char soilwaterf[255];   /* soil water file, XF, 24 June, 2003*/
char latitudef[255];    /* Latitude data file */
char longitudef[255];   /* Longitude data file */
char tmaxf[255];        /* Temperature Maxmium daily data file */
char tminf[255];        /* Temperature Minmium daily data file */
char precipf[255];      /* precipitation daily file */
char rehumf[255];       /* Relative humidity daily file */
char sradf[255];        /* soliar radiation daily file */
char leafmaxf[255];     //zfm3.26


int long cycle_display;

int long lin_offset;    /* Line offset (pci coord. system) */

int long pix_offset;    /* Pixel offset (pci coord. system) */
 
float factor;
char outf_name[255];

char output_NPP;
char output_GPP;
char output_res;
char output_eva;
char output_tra;
char output_site;

