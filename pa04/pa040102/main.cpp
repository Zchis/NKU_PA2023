#include<pthread.h>
#include<iostream>
#include<malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include<random>
#include <sys/time.h>
#include<immintrin.h>
using namespace std;
const int N=1000;
const int n=N;
float A[N][N];
void m_reset()
{
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<i;j++)
			A[i][j]=0;
		A[i][i]=1.0;
		for(int j=i+1;j<N;j++)
			A[i][j]=rand();
	}
	for(int k=0;k<N;k++)
		for(int i=k+1;i<N;i++)
			for(int j=0;j<N;j++)
				A[i][j]+=A[k][j];
}
typedef struct{
    int k;
    int t_id;
}threadParam_t;
void *threadFunc(void *param){
    threadParam_t *p=(threadParam_t*)param;
    int k=p->k;
    int t_id=p->t_id;
    int i=k+t_id+1;

    __m128 vaik = _mm_set1_ps(A[i][k]);
    int j = k + 1;
    for(;j+4<=N;j+=4){
        __m128 vakj = _mm_loadu_ps(&A[k][j]);
        __m128 vaij = _mm_loadu_ps(&A[i][j]);
        __m128 vx = _mm_set1_ps(0);
        vx = _mm_mul_ps(vakj, vaik);//乘法
        vaij = _mm_sub_ps(vaij, vx);//减法
        _mm_store1_ps(&A[i][j], vaij);//加载回内存
    }
    for (; j < N; j++) {
        A[i][j] = A[i][j] - A[k][j] * A[i][k];
    }

    A[i][k]=0;
    pthread_exit(NULL);
}
int main(){
    m_reset();
    struct timeval my_start;
	struct timeval my_end;//clock
	float timecount;
	gettimeofday(&my_start,NULL);
    for(int k=0;k<N;k++){
/*
        for(int j=k;j<N;j++){
            A[k][j]=A[k][j]/A[k][k];
        }
*/
        int j=k+1;
        __m128 vt = _mm_set1_ps(A[k][k]);

        for(;j+4<=N;j+=4){
            __m128 va = _mm_loadu_ps(&A[k][j]);
            va = _mm_div_ps(va, vt);
            _mm_store1_ps(&A[k][j], va);
        }
        for(;j<N;j++){
            A[k][j]=A[k][j]/A[k][k];
        }

        A[k][k]=1.0;

        int worker_count=N-1-k;
        pthread_t* handles=(pthread_t*)malloc(worker_count*sizeof(pthread_t));
        threadParam_t* param=(threadParam_t*)malloc(worker_count*sizeof(threadParam_t));


        for(int t_id=0;t_id<worker_count;t_id++){
            param[t_id].k=k;
            param[t_id].t_id=t_id;
        }

        for(int t_id=0;t_id<worker_count;t_id++){
            pthread_create(&handles[t_id],NULL,threadFunc,&param[t_id]);
        }

        for(int t_id=0;t_id<worker_count;t_id++){
            pthread_join(handles[t_id],NULL);
        }
    }
	gettimeofday(&my_end,NULL);
    timecount+=(my_end.tv_sec-my_start.tv_sec)*1000000+my_end.tv_usec-my_start.tv_usec;
	cout<<timecount/1000<<endl;
    return 0;
}
