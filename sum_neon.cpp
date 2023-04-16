#include<iostream>
#include<arm_neon.h>
#include<sys/timeb.h>
using namespace std;
float sum_ord(float *array,int len){
    float sum=0.0;
    for(int i=0;i<len;i++){
        sum+=array[i];
    }
    return sum;
}
float sum_neon(float* array,int len){
    float sum=0.0;
    float32x4_t sum_vec = vdupq_n_f32(0.0);
    for(int i=len-4;i>=0;i-=4){
        float32x4_t temp=vld1q_f32(array+i);
        sum_vec = vaddq_f32(sum_vec,temp);
    }
    float sum=vgetq_lane_f32(sum_vec,1)+vgetq_lane_f32(sum_vec,2)+vgetq_lane_f32(sum_vec,3)+vgetq_lane_f32(sum_vec,0);
    for(int i=0;i%4!=(len%4);i++){
        sum += array[i];
    }
    return sum;
}
int main(){
    int N;
    cin>>N;
    float* a=new float[N];
    for(int i=0;i<N;i++){
        a[i]=i;
    }
    struct  timeval   tv_begin,tv_end;
    float time=0.0;

    gettimeofday(tv_begin,NULL);
    sum_ord(a,N);
    gettimeofday(tv_end,NULL);
    time = (tv_end.tv_sec-tv_begin.tv_sec)*1000000+(tv_end.tv_usec-tv_end.tv_usec);
    cout<<"the time of ord is "<<time<<" ms "<<endl;

    gettimeofday(tv_begin,NULL);
    sum_neon(a,N);
    gettimeofday(tv_end,NULL);
    time = (tv_end.tv_sec-tv_begin.tv_sec)*1000000+(tv_end.tv_usec-tv_end.tv_usec);
    cout<<"the time of ord is "<<time<<" ms "<<endl;

    return 0;
}