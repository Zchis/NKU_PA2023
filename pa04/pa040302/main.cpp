#include<pthread.h>
#include<iostream>
#include<malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include<random>
#include <sys/time.h>
#include <semaphore.h>
#include<immintrin.h>
using namespace std;
const int N=1000;
const int n=N;
float A[N][N];
const double eps = 1e-6;
const int NUM_THREADS=7;
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
    int t_id;
}threadParam_t;

sem_t sem_leader;
sem_t sem_Division[NUM_THREADS-1];
sem_t sem_Elimination[NUM_THREADS-1];

void *threadFunc(void *param){
    threadParam_t *p=(threadParam_t*)param;
    int t_id=p->t_id;

    for(int k=0;k<n;++k){
        if(t_id==0){
            /*
            for(int j=k+1;j<n;j++){
                A[k][j]=A[k][j]/A[k][k];
            }
            A[k][k]=1.0;
            */
            int j=k+1;
            __m128 vt = _mm_set1_ps(A[k][k]);

            for(;j+4<N;){
                __m128 va = _mm_loadu_ps(&A[k][j]);
                va = _mm_div_ps(va, vt);
                _mm_store1_ps(&A[k][j], va);
                j=j+4;
            }

            for(;j<N;j++){
                A[k][j]=A[k][j]/A[k][k];
            }
            A[k][k] = 1.0;
        }
        else{
            sem_wait(&sem_Division[t_id-1]);
        }
        if(t_id==0){
            for(int i=0;i<NUM_THREADS-1;i++){
                sem_post(&sem_Division[i]);
            }
        }
        for(int i=k+1+t_id;i<n;i+=NUM_THREADS){
            /*
            for(int j=k+1;j<n;j++){
                A[i][j]-=A[i][k]*A[k][j];
            }
            A[i][k]=0.0;
            */
            __m128 vaik = _mm_set1_ps(A[i][k]);
            int j = k + 1;
            for(;j+4<=N;j+=4){
                __m128 vakj = _mm_loadu_ps(&A[k][j]);
                __m128 vaij = _mm_loadu_ps(&A[i][j]);
                __m128 vx = _mm_set1_ps(0);
                vx = _mm_mul_ps(vakj, vaik);//�˷�
                vaij = _mm_sub_ps(vaij, vx);//����
                _mm_store1_ps(&A[i][j], vaij);//���ػ��ڴ�
            }
            for (; j < N; j++) {
                A[i][j] = A[i][j] - A[k][j] * A[i][k];
            }
            A[i][k]=0.0;
        }
        if(t_id==0){
            for(int i=0;i<NUM_THREADS-1;i++){
                sem_wait(&sem_leader);
            }
            for(int i=0;i<NUM_THREADS-1;i++){
                sem_post(&sem_Elimination[i]);
            }

        }
        else{
            sem_post(&sem_leader);
            sem_wait(&sem_Elimination[t_id-1]);
        }
    }
    pthread_exit(NULL);
}


int main()
{
    m_reset();

    struct timeval my_start;
	struct timeval my_end;//clock
	float timecount;
	gettimeofday(&my_start,NULL);

    sem_init(&sem_leader,0,0);
    for(int i=0;i<NUM_THREADS-1;i++){
        sem_init(&sem_Division[i],0,0);
        sem_init(&sem_Elimination[i],0,0);
    }

    pthread_t handles[NUM_THREADS];
    threadParam_t param[NUM_THREADS];
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        param[t_id].t_id=t_id;
        pthread_create(&handles[t_id],NULL,threadFunc,&param[t_id]);
    }

    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        pthread_join(handles[t_id],NULL);
    }
    sem_destroy(&sem_leader);
        for(int i=0;i<NUM_THREADS-1;i++){
        sem_destroy(&sem_Division[i]);
        sem_destroy(&sem_Elimination[i]);
    }


    gettimeofday(&my_end,NULL);
    timecount+=(my_end.tv_sec-my_start.tv_sec)*1000000+my_end.tv_usec-my_start.tv_usec;
	cout<<timecount/1000<<endl;

    return 0;
}
