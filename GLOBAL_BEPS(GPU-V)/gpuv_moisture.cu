
#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 


/* This subroutine is to estimate the vertical redistribution of soil water*/
__device__ void gpuv_moisture(int long pix, float depth3, float B[], float K[], float Psi[],
                              struct hy_c1 HY1[], struct hy_c2 HY2[], struct hy_c3 HY3[], float *theta_o, float *theta_n1,float *theta_n2, float *theta_n3)
{

float X[10][10],Y[10];
float x1[10],x2[10],x3[10];
//float Flux[5];// downward water across layer boundaries 
float m[5];
float n[5];

float bperm;         // =0 for impermeable (no drainage)  lower bc; =1 for free drainage lower bc
float hsoim[5]; //vertical distance between cernters of layers 
float wsoim[5]; // interpolated moisture values at layer boundaries 
float wsoia[5]; // interpolated moisture values at layer boundaries  
float wsoib[5]; // interpolated moisture values at layer boundaries  

float weim[5]; // interpolated moisture values at layer boundaries  
float weip[5]; // interpolated moisture values at layer boundaries  

int i;  //intermideate variable
float e[5],f[5],g[5];
float a[5],b[5],zbex;
float bwn1[5],bwn[5];
float rimp ;// the implicit f
//float rhs[10];
float dt;

//float rhow=1000.0; //the density of water
int km1;

int layer1;
float surface_plus;
float kk=1;
float d_n,d_o;
float SW[5], SS[5],depth[5],SF[5];//
float interval;
float thetan[5];


bperm=0;
interval=1.0;
layer1=layer;
rimp=1.0;//1: complete implicit solution, 0 explicit solution

SW[1]=HY1[pix].WP; SF[1]=HY1[pix].FC;  SS[1]=HY1[pix].PR;
SW[2]=HY2[pix].WP; SF[2]=HY2[pix].FC;  SS[2]=HY2[pix].PR;
SW[3]=HY3[pix].WP; SF[3]=HY3[pix].FC;  SS[3]=HY3[pix].PR;
depth[1]=depth1;depth[2]=depth2; depth[3]=depth3;


  hsoim[1]=0.5*depth[1];   weim[1]=0;   weip[1]=1.0;  wsoim[1]=theta_o[1];
  wsoia[1]=__min(wsoim[1],1.0); 
  wsoib[1]=__min(wsoim[1],1.0); 

  for(i=2;i<=layer1;i++){
  hsoim[i]=0.5*(depth[i-1]+depth[i]);
  //weim[i]=0.5*depth[i-1]/hsoim[i];    //ORGIS:   weim[i]=0.5*depth[i]/hsoim[i];  
 
  weim[i]=0.5*depth[i]/hsoim[i];    //2018-03-02恢复 
  weip[i]=1-weim[i];




  wsoim[i]=weim[i]*theta_o[i-1]+weip[i]* theta_o[i];
  wsoia[i]=__min(wsoim[i],1.0);
  wsoib[i]=__min(wsoim[i],1.0);
  }  // the end of i loop
  
  i=layer1+1;
  
  hsoim[i]=0.5 *depth[i-1];
  weim[i]=1.0;
  weip[1]=0;
   
  wsoim[i]=theta_o[i-1];
  wsoia[i]=__min(wsoim[i],1.0);
  wsoib[i]=__min(wsoim[i],1.0);

   
 e[1]=0; f[1]=0; g[1]=0;
 	  
 for (i=2;i<=layer1;i++) {
 a[i]=weim[i]*K[i-1]+weip[i]*K[i];
 b[i]=weim[i]*K[i-1]*Psi[i-1]*B[i-1]+weip[i]*K[i]*Psi[i]*B[i];
 zbex=weim[i]*B[i-1]+weip[i]*B[i];   // B 
 m[i]=2*zbex+3;  
 n[i]=zbex+2; 

 bwn1[i]=b[i]*pow(wsoib[i],(n[i]-1));

 bwn[i]=bwn1[i]*wsoib[i];



 e[i]=a[i]*(-1+rimp*m[i])*pow(wsoia[i],m[i])
	 +((1-rimp)*bwn[i]-rimp*n[i]*bwn1[i]*wsoib[i])*(theta_o[i]-theta_o[i-1])/hsoim[i];
 
 f[i]=-rimp*a[i]*m[i]*pow(wsoia[i],(m[i]-1))
	  +rimp*n[i]*bwn1[i]*(theta_o[i]-theta_o[i-1])/hsoim[i];
 



 /*
 e[i]=a[i]*(1-rimp*m[i])*pow(wsoia[i],m[i])
	 +((1-rimp)*bwn[i]-rimp*n[i]*bwn1[i]*wsoib[i])*(theta_o[i]-theta_o[i-1])/hsoim[i];

 f[i]=rimp*a[i]*m[i]*pow(wsoia[i],(m[i]-1))
	 +rimp*n[i]*bwn1[i]*(theta_o[i]-theta_o[i-1])/hsoim[i];


*/

 g[i]= rimp*bwn[i]; 
	
 } // the end of i loop

i=layer1+1;   // for the deepest layer
a[i]= K[i-1];
b[i]= K[i-1]*Psi[i-1]*B[i-1];
m[i]= 2*B[i-1]+3;
n[i]= B[i-1]+2;
d_o=rimp*a[i]*m[i]*pow(wsoia[i],m[i]-1)*theta_o[i-1]*bperm;
d_n=rimp*a[i]*m[i]*pow(wsoia[i],m[i]-1)*bperm;

e[i]=0;//-(a[i]*pow(wsoia[i],m[i])*bperm-d_o*kk); 
 //the considertion of vertical drainage from the bottom of soil profile
f[i]=0;
g[i]=0;



/*********************************************************************************************************/
for (i=1;i<=layer1;i++){
dt=interval/(SS[i]*depth[i]);
X[i][1]=dt*(f[i]*0.5*depth[i]/hsoim[i]-g[i]/hsoim[i]);
Y[i]=theta_o[i]+dt*(e[i+1]-e[i]); 



if(i<layer1){
km1=__max(i-1,1);
X[i][2]=1+dt*(-f[i+1]*0.5*depth[i+1]/hsoim[i+1]+f[i]*0.5*depth[km1]/hsoim[i]+g[i+1]/hsoim[i+1]+g[i]/hsoim[i]);
X[i][3]=  dt*(-f[i+1]*0.5*depth[i]/hsoim[i+1]-g[i+1]/hsoim[i+1]);
}  // end of i<layer

else if (i==layer){     //for the deepest layer 
dt=interval/(SS[i]*depth[i]);
X[i][2]=1.0+dt*(-f[i+1]+f[i]*0.5*depth[i-1]/hsoim[i]+g[i]/hsoim[i]);
X[i][2]=X[i][2]+dt*d_n*kk;           //*rimp*m[i]*a[i+1]*pow(wsoia[i+1],m[i]-1)*bperm;
//

X[i][3]=0;
}  // the end of else
}  //the end  of i loop

for(i=0;i<10;i++){
x1[i]=X[1][i];
x2[i]=X[2][i];
x3[i]=X[3][i];
}



 
 gpusolution11(thetan,x1,x2,x3,Y);

/*
for(i=layer1;i>=1;i--){
if(thetan[i]>1.0) {
if(i>1){
thetan[i-1]=thetan[i-1]+(thetan[i]-1)*SS[i]*depth[i]/(SS[i-1]*depth[i-1]);
thetan[i]=1.0;
	}
else {
//surface_plus=(theta_n[i]-1.0)*SS[i]*depth[i];	
//theta_n[i]=1.0;	
}

}   // the end of theta_n[i]>1.0

//if(theta_n[i]<(SW[i]/SS[i]*0.75))theta_n[i]=SW[i]/SS[i]*0.75;  
//系数0.5为2012-10-15日加    //2017-06-02发现上面的该限制有问题，会导致模拟的土壤水分偏高， 取消该限制  

}  // the end of i loop




if(thetan[1]>1.0){
surface_plus=(thetan[1]-1.0)*SS[1]*depth[1];
thetan[1]=1.0;
}
else  surface_plus=0;

*/


*theta_n1=thetan[1];

*theta_n2=thetan[2];

*theta_n3=thetan[3];

}  // the end of subroutine