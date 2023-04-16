#include<iostream>
#include<arm_neon.h>
#include<random>

#include<stdio.h>
#include<sys/time.h>


using namespace std;
int main(){
  struct timespec sts,ets;
  const int n=1000;
  float A[n][n];
  for(int i=0;i<n;i++){
    A[i][i]=1.0;
    for(int j=0;j<n;j++){
      A[i][j]=rand();
    }
  }
  for(int k=0;k<n;k++){
    for(int i=k+1;i<n;i++){
      for(int j=0;j<n;j++){
        A[i][j]+=A[k][j];
      }
    }
  }
timespec_get(&sts, TIME_UTC);

  for(int k=0;k<n;k++){
    float32x4_t vt = vmovq_n_f32(A[k][k]);
    for(int j=k+1;j<n;j+=4){
      float32x4_t va = vld1q_f32(&A[k][j]);
      va = vdivq_f32(va, vt);
      vst1q_f32(&A[k][j], va);
    }

    for(int i=k+1;i<n;i++){
      float32x4_t vaik=vmovq_n_f32(A[i][k]);
      int j=k+1;
      for(;j+4<=n;j+=4){
        float32x4_t vakj = vld1q_f32(&A[k][j]);
        float32x4_t vaij = vld1q_f32(&A[i][j]);
        float32x4_t vx = vmovq_n_f32(0);
        vx = vmulq_f32(vakj, vaik);
        vaij = vsubq_f32(vaij, vx);
        vst1q_f32(&A[i][j], vaij);
      }
      for(;j<n;j++){
        A[i][j]=A[i][j]-A[k][j]*A[i][k];
      }
    }
  }
timespec_get(&ets, TIME_UTC);
time_t dsec=ets.tv_sec-sts.tv_sec;
long dnsec=ets.tv_nsec-sts.tv_nsec;

if(dnsec<0){
dsec--;
dnsec+=1000000000ll;
}
printf("%ld.%09lds\n",dsec,dnsec);
cout<<"Over!"<<endl;
}
