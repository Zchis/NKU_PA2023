#include<iostream>
#include<arm_neon.h>
#include<sys/timeb.h>
using namespace std;
int main()
{
    struct  timeval   tv_begin,tv_end;
    double time;
    int N;
    cin>>N;
    float32_t** A=new float32_t*[N];
    for(int i=0;i<N;i++){
        A[i]=new float32_t[N];
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            A[i][j]=i+j+1;
        }
    }
    
    gettimeofday(tv_begin,NULL);
    for(int k=0;k<N;k++){
        for(int j=k;j<N;j++){
            A[k][j]=A[k][j]/A[k][k];
        }
        A[k][k]=1.0;
        for(int i=k+1;i<N;i++){
            for(int j=k+1;j<N;j++){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
        }
    }

    gettimeofday(tv_end,NULL);
    time = (tv_end.tv_sec-tv_begin.tv_sec)*1000000+(tv_end.tv_usec-tv_end.tv_usec);
    cout<<"the time of a way is "<<time<<" ms "<<endl;

    float32x4_t a,b,c,d;
    gettimeofday(tv_begin,NULL);
    for(int k=0;k<N;k++){
        a=vdup_n_f32(A[k][k]);
        for(int j=N-4;j>=k;j-=4){
            b=vld1q_f32(A[k]+j);
            c=
            A[k][j]=A[k][j]/A[k][k];
        }
        A[k][k]=1.0;
        for(int i=k+1;i<N;i++){
            for(int j=k+1;j<N;j++){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
        }
    }

    gettimeofday(tv_end,NULL);
    time = (tv_end.tv_sec-tv_begin.tv_sec)*1000000+(tv_end.tv_usec-tv_end.tv_usec);
    cout<<"the time of a way is "<<time<<" ms "<<endl;

    return 0;
}