#include<pthread.h>
#include<iostream>
#include<malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include<random>
#include <sys/time.h>
using namespace std;
const int N=1000;
float A[N][N];
const double eps = 1e-6;
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

    for(int j=k+1;j<N;j++){
        A[i][j]=A[i][j]-A[i][k]*A[k][j];
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
        for(int j=k;j<N;j++){
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
