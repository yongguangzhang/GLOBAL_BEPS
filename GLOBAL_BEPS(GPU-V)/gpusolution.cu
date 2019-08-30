
#ifndef   MY_H_FILE      
#define   MY_H_FILE       
#include "gpubeps.cuh"
#endif 

__device__ void gpusolution11(float xx[],float x1[],float x2[],float x3[],float Y[])
{
int i,j;
float alpha[10],gamma[10],a[10],b[10],c[10];
int num;
num=layer;

float X[10][10];

for(i=0;i<10;i++){
	X[1][i]=x1[i];
	X[2][i]=x2[i];
	X[3][i]=x3[i];
}





for(i=1;i<=num;i++){
a[i]=X[i][1];b[i]=X[i][2];c[i]=X[i][3];
}
 
alpha[1]=1.0/b[1];
gamma[1]=c[1]*alpha[1];

for(i=2;i<num;i++){
 
alpha[i]=1/(b[i]-a[i]*gamma[i-1]);
gamma[i]=c[i]*alpha[i];
}


xx[1]=Y[1]*alpha[1];

for(i=2;i<num;i++) xx[i]=(Y[i]-a[i]*xx[i-1])*alpha[i];

xx[num]=(Y[num]-a[num]*xx[num-1])/(b[num]-a[num]*gamma[num-1]);
 
for(i=1;i<num;i++){
j=num-i;
xx[j]=xx[j]-gamma[j]*xx[j+1];
}

}
