#include <iostream>
#include<windows.h>
#include<xmmintrin.h>  //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX¡¢AVX2
using namespace std;

int main()
{
    int N;
    cin>>N;
    float** A=new float*[N];
    for(int i=0;i<N;i++){
        A[i]=new float[N];
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            A[i][j]=i+j+1;
        }
    }
    double time;
    LARGE_INTEGER freq_;
    QueryPerformanceFrequency(&freq_);

    LARGE_INTEGER begin_time;
    LARGE_INTEGER end_time;

    //串行算法
    /*QueryPerformanceCounter(&begin_time);
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
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart-begin_time.QuadPart)/(double)freq_.QuadPart;
    cout<<"the N is "<<N<<endl;
    cout<<"the time of ord costs "<<time*1000<<" ms"<<endl;*/

    //SSE指令优化
    __m128 a,b,c,d;
    /*QueryPerformanceFrequency(&freq_);
    QueryPerformanceCounter(&begin_time);
    for(int k=0;k<N;k++){
        float temp[4]={A[k][k],A[k][k],A[k][k],A[k][k]};
        a=_mm_loadu_ps(temp);
        for(int j=N-4;j>=k;j-=4){
            b=_mm_loadu_ps(A[k]+j);
            c=_mm_div_ps(b,a);
            _mm_store_ps(A[k]+j,c);
        }
        if(k%4!=(N%4)){
            for(int j=k;j%4!=(N%4);j++){
                A[k][j]=A[k][j]/A[k][k];
            }
        }
        A[k][k]=1.0;
        for(int i=k+1;i<N;i++){
            for(int j=k+1;j<N;j++){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart-begin_time.QuadPart)/(double)freq_.QuadPart;
    cout<<"the N is "<<N<<endl;
    cout<<"1.the time of SSE costs "<<time*1000<<" ms"<<endl;*/

    /*QueryPerformanceFrequency(&freq_);
    QueryPerformanceCounter(&begin_time);
    for(int k=0;k<N;k++){
        for(int j=k+1;j<N;j++){
            A[k][j]=A[k][j]/A[k][k];
        }
        A[k][k]=1.0;
        for(int i=k+1;i<N;i++){
        float temp1[4]={A[i][k],A[i][k],A[i][k],A[i][k]};
        a=_mm_loadu_ps(temp1);
            for(int j=N-4;j>k;j-=4){
                b=_mm_loadu_ps(A[i]+j);
                c=_mm_loadu_ps(A[i]+j);
                d=_mm_sub_ps(b,_mm_mul_ps(a,c));
                _mm_store_ps(A[i]+j,d);
            }
            for(int j=k+1;j%4!=(N%4);j++){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart-begin_time.QuadPart)/(double)freq_.QuadPart;
    cout<<"the N is "<<N<<endl;
    cout<<"2.the time of SSE costs "<<time*1000<<" ms"<<endl;*/


    //SSE不对齐
    /*QueryPerformanceFrequency(&freq_);
    QueryPerformanceCounter(&begin_time);
    for(int k=0;k<N;k++){
        float temp[4]={A[k][k],A[k][k],A[k][k],A[k][k]};
        a=_mm_loadu_ps(temp);
        for(int j=N-4;j>=k;j-=4){
            b=_mm_loadu_ps(A[k]+j);
            c=_mm_div_ps(b,a);
            _mm_store_ps(A[k]+j,c);
        }
        if(k%4!=(N%4)){
            for(int j=k;j%4!=(N%4);j++){
                A[k][j]=A[k][j]/A[k][k];
            }
        }
        A[k][k]=1.0;
        for(int i=k+1;i<N;i++){
        float temp1[4]={A[i][k],A[i][k],A[i][k],A[i][k]};
        a=_mm_loadu_ps(temp1);
            for(int j=N-4;j>k;j-=4){
                b=_mm_loadu_ps(A[i]+j);
                c=_mm_loadu_ps(A[i]+j);
                d=_mm_sub_ps(b,_mm_mul_ps(a,c));
                _mm_store_ps(A[i]+j,d);
            }
            for(int j=k+1;j%4!=(N%4);j++){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart-begin_time.QuadPart)/(double)freq_.QuadPart;
    cout<<"the N is "<<N<<endl;
    cout<<"the time of SSE costs "<<time*1000<<" ms"<<endl;*/

    //SSE对齐
    QueryPerformanceFrequency(&freq_);
    QueryPerformanceCounter(&begin_time);
    for(int k=0;k<N;k++){
        float temp[4]={A[k][k],A[k][k],A[k][k],A[k][k]};
        a=_mm_loadu_ps(temp);
        for(int j=N-4;j>=k;j-=4){
            b=_mm_load_ps(A[k]+j);
            c=_mm_div_ps(b,a);
            _mm_store_ps(A[k]+j,c);
        }
        if(k%4!=(N%4)){
            for(int j=k;j%4!=(N%4);j++){
                A[k][j]=A[k][j]/A[k][k];
            }
        }
        A[k][k]=1.0;
        for(int i=k+1;i<N;i++){
        float temp1[4]={A[i][k],A[i][k],A[i][k],A[i][k]};
        a=_mm_loadu_ps(temp1);
            for(int j=N-4;j>k;j-=4){
                b=_mm_load_ps(A[i]+j);
                c=_mm_load_ps(A[i]+j);
                d=_mm_sub_ps(b,_mm_mul_ps(a,c));
                _mm_store_ps(A[i]+j,d);
            }
            for(int j=k+1;j%4!=(N%4);j++){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart-begin_time.QuadPart)/(double)freq_.QuadPart;
    cout<<"the N is "<<N<<endl;
    cout<<"the time of SSE costs "<<time*1000<<" ms"<<endl;


    //AVX指令优化
    /*__m256 a1,b1,c1,d1;
    QueryPerformanceFrequency(&freq_);
    QueryPerformanceCounter(&begin_time);
    for(int k=0;k<N;k++){
        float temp[8]={A[k][k],A[k][k],A[k][k],A[k][k],A[k][k],A[k][k],A[k][k],A[k][k]};
        a1=_mm256_loadu_ps(temp);
        for(int j=N-8;j>=k;j-=8){
            b1=_mm256_loadu_ps(A[k]+j);
            c1=_mm256_div_ps(b1,a1);
            _mm256_store_ps(A[k]+j,c1);
        }
        if(k%8!=(N%8)){
            for(int j=k;j%8!=(N%8);j++){
                A[k][j]=A[k][j]/A[k][k];
            }
        }
        A[k][k]=1.0;
        for(int i=k+1;i<N;i++){
            float temp1[8]={A[i][k],A[i][k],A[i][k],A[i][k],A[i][k],A[i][k],A[i][k],A[i][k]};
            a1=_mm256_loadu_ps(temp1);
            for(int j=N-8;j>k;j-=8){
                b1=_mm256_loadu_ps(A[i]+j);
                c1=_mm256_loadu_ps(A[i]+j);
                d1=_mm256_sub_ps(b1,_mm256_mul_ps(a1,c1));
                _mm256_store_ps(A[i]+j,d1);
            }
            for(int j=k+1;j%8!=(N%8);j++){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart-begin_time.QuadPart)/(double)freq_.QuadPart;
    cout<<"the N is "<<N<<endl;
    cout<<"the time of AVXu costs "<<time*1000<<" ms"<<endl;*/

    //AVX指令优化，对齐
    /*QueryPerformanceFrequency(&freq_);
    QueryPerformanceCounter(&begin_time);
    for(int k=0;k<N;k++){
        float temp[8]={A[k][k],A[k][k],A[k][k],A[k][k],A[k][k],A[k][k],A[k][k],A[k][k]};
        a1=_mm256_load_ps(temp);
        for(int j=N-8;j>=k;j-=8){
            b1=_mm256_load_ps(A[k]+j);
            c1=_mm256_div_ps(b1,a1);
            _mm256_store_ps(A[k]+j,c1);
        }
        if(k%8!=(N%8)){
            for(int j=k;j%8!=(N%8);j++){
                A[k][j]=A[k][j]/A[k][k];
            }
        }
        A[k][k]=1.0;
        for(int i=k+1;i<N;i++){
        float temp1[8]={A[i][k],A[i][k],A[i][k],A[i][k],A[i][k],A[i][k],A[i][k],A[i][k]};
        a1=_mm256_loadu_ps(temp1);
            for(int j=N-8;j>k;j-=8){
                b1=_mm256_load_ps(A[i]+j);
                c1=_mm256_load_ps(A[i]+j);
                d1=_mm256_sub_ps(b1,_mm256_mul_ps(a1,c1));
                _mm256_store_ps(A[i]+j,d1);
            }
            for(int j=k+1;j%8!=(N%8);j++){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart-begin_time.QuadPart)/(double)freq_.QuadPart;
    cout<<"the N is "<<N<<endl;
    cout<<"the time of AVX costs "<<time*1000<<" ms"<<endl;*/



    //循环展开
    /*QueryPerformanceCounter(&begin_time);
    for(int k=0;k<N;k++){
        for(int j=N-4;j>=k;j-=4){
            A[k][j]=A[k][j]/A[k][k];
            A[k][j+1]=A[k][j+1]/A[k][k];
            A[k][j+2]=A[k][j+2]/A[k][k];
            A[k][j+3]=A[k][j+3]/A[k][k];
        }
        if(k%4!=(N%4)){
            for(int j=k;j%4!=(N%4);j++){
                A[k][j]=A[k][j]/A[k][k];
            }
        }
        A[k][k]=1.0;
        for(int i=k+1;i<N;i++){
            for(int j=N-4;j>k;j-=4){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
                A[i][j+1]=A[i][j+1]-A[i][k]*A[k][j+1];
                A[i][j+2]=A[i][j+2]-A[i][k]*A[k][j+2];
                A[i][j+3]=A[i][j+3]-A[i][k]*A[k][j+3];
            }
            for(int j=k+1;j%4!=(N%4);j++){
                A[i][j]=A[i][j]-A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart-begin_time.QuadPart)/(double)freq_.QuadPart;
    cout<<"the N is "<<N<<endl;
    cout<<"the time of loopunroll costs "<<time*1000<<" ms"<<endl;*/

    return 0;
}
