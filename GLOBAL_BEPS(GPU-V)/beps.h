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
/* Define Constants */

#define NOERROR		 0
#define ERROR		 1
#define PI		     3.1415926


//#define nlines       2090         //计算区域的行数 
//#define npixels      4950    //计算区域的列数





#define nlines       5         //计算区域的行数 
#define npixels      2069100    //计算区域的列数


/*

#define nlines       209         //计算区域的行数 
#define npixels      49500    //计算区域的列数

#define nlines       418                    //计算区域的行数 
#define npixels      24750                     //计算区域的列数


#define nlines       5         //计算区域的行数 
#define npixels      2069100    //计算区域的列数

#define nlines       10       //计算区域的行数 
#define npixels      1034550    //计算区域的列数


 */

#define layer        3
#define jday_start   1                // Julian day start 
#define jday_end     365              //lixuansong


#define SIZEB		 61
#define SIZEG 	    61
#define SIZEX 		60#define SIZEZ 		 70  

#define RTIMES 		 1
#define RPERIOD 	 3600	// in sec 


//bool InitCUDA();
void solution11(float xx[]);
void v_moisture(int long pix,  
				float B[],
				float K[],
				float Psi[],
struct hy_c1 HY1[],
struct hy_c2 HY2[],
struct hy_c3 HY3[],
	float *theta_o,
	float *theta_n);
void inter(
		   float  b[],float g[],float x[],float z[]);
void melt (
		   float b[],float g[],float x[],float z[]);
void rs(
		int long pix,
		float b[],float g[],float x[],float z[],float rr[],
struct hy_c1 HY1[],
struct hy_c2 HY2[], 
struct hy_c3 HY3[]);
float penmon(float tair,float vpd,float wrad,float acond,float scond);
void readxx(int long pix,float x[],struct xvalue xx[]);
void writexx(int long pix,float x[],struct xvalue xx[]);

void doflux(
			int long pix,
			float b[],float g[],float x[],float z[],
			float co,
			float soilw1,float soilw2,float soilw3,

struct hy_c1 HY1[],
struct hy_c2 HY2[],
struct hy_c3 HY3[],
	float *soilwn1,float*soilwn2,float*soilwn3);
void readb_init(float b[]); 
void zeroxx1(int long pix,float x[],struct xvalue xx[]);
void zeroxx(int long pix,float x[],struct xvalue xx[]);
void readlc(int long line,int short *lc,int short *rcode);
void readsoildata(long line,struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[],int short *rcode);
void readlaimax(int long line,struct xvalue xx[],int short *lc,int short *rcode,int tmpyear);


void lw(float FC1,float FC2, float FC3, float sw1, float  sw2, float sw3, float *lw1,float *lw2, float *lw3);


//void readNdep(int long line,float Ndep[],int short *rcode,int tmpyear);

void readclim(long year,long jday,long line,struct climatedata sdat[],int short *rcode);
void readlai(int long line,float *lai,int short jday,int short *rcode,long tmpyear);
void readb(float *b, short int lc_p);
void coef_c(long pix,short lc,float coef[],struct hy_c2 HY2[]);
void model(long jday,int long pix,float lat_p,float lai_p,  short int lc_p,float sw_o1,float sw_o2,float sw_o3,struct hy_c1 HY1[],struct hy_c2 HY2[],struct hy_c3 HY3[],
float     x[],
float     b[],
struct     climatedata sdat[],
struct     xvalue xx[],
float      co);
void preInitials(float lambd,float lambddd[],float  NPP,float coef1[],float initialValues[]);

void soilresp(
			  float coef[],
float  lambda,float Cstem1,float Cleaf1,float Ccroot1,float Cfroot1,float Ccd1,float Csmd1,float Cssd1,float Cfmd1 ,float Cfsd1,float Csm1 ,float Cm1 ,float Cs1 ,float Cp1,
float  initialValues[]);

void setx(short int lc_p, float awc_p, float x[], float b[]);
void zcomp(
		   long jday,
		   int long pix,
		   float lat_p,
		   short  int lc_p,
		   float awc_p,
		   float lai_p,
		   float b[],float   x[],float z[],
struct climatedata sdat[]);
int rad_ssl(
			float sg,			        
			float cos_theta,
			float lai_p,	  
			float omega,		  
			float *ssun,		  
			float *sshade);
int net_shortwave(
				  float sg,			     
				  float cos_theta,
				  float lai_p,		     
				  float omega,		     
				  float *ssun,		    
				  float *sshade,
				  float lai_u,		    
				  float *n_ssun,
				  float *n_sshade,
				  float *n_sunder,		 
				  float *n_sground);
int net_longwave(
				 float lai_o,
				 float lai_u,
				 float omega,
				 float es,
				 float ta,
				 float *m_lnet_o,
				 float *m_lnet_u,		             
				 float *lnet_g);		
void carbon(  
			float b[],float g[],float x[],float z[],float co);
int farq_psn(float b[],float z[],float g[],float co);
__device__ float Max(float a,float b);
__device__ float Min(float a,float b);

/*
void readconf();
void readb_init();
void pixtolon();
void readlc();
void readlatit();
void readlongi();
void readsoildata();
void readsnow();
void readlaimax	();
void readlai();
void readclim();
void readb();
void model();
void coef_c();
void preInitials();
void soilresp();
void v_moisture();
void readsoildata();
void setx();
void zeroxx1();   //zfm4.21
void zeroxx();
void readxx();
void writexx();

void zcomp();
void doflux();
void lw();

void inter();
void melt();
void ra();
void rs();
void rsoil();
void carbon();

int net_longwave();
int net_shortwave();
int farq_psn();
int rad_ssl();

float penman();

void display();
void displaynpp();
void solution11();

*/

// Declare structures 

struct climatedata
	{
	int long pix;
    float tmax;
	float tmin;
	float rehum;
	float srad;
	float precip;
  };

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

struct hy_c3
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
	float x10;		/* root C */
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


	float x31;    
	float x32;    
	float x33;    
	float x34;    


	float x41;    
	float x42;    
	float x43;    
	float x44;    

};


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

//lixuansong:全局变量转局部变量
//float *z;
//float depth[10];
float X[10][10],Y[10];
//float baita;

//float soilw1_new, soilw2_new, soilw3_new;


//gpu
__device__ void gpusolution11(float xx[],float x1[],float x2[], float x3[],float Y[]);