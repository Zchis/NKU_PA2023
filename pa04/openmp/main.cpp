#include <iostream>
#include <random>
#include<omp.h>
#include<windows.h>
#include<xmmintrin.h>  //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX¡¢AVX2
#define THREADS_NUM 4

using namespace std;

float** A;
int N;
void init() {
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
			A[i][j] = 0.0;
		A[i][i] = 1.0;
		for (int j = i + 1; j < N; j++)
			A[i][j] = rand();
	}
	for (int k = 0; k < N; k++)
		for (int i = k + 1; i < N; i++)
			for (int j = 0; j < N; j++)
				A[i][j] += A[k][j];
}


LARGE_INTEGER freq_;
LARGE_INTEGER begin_time;
LARGE_INTEGER end_time;
int main()
{
    //omp_set_num_threads(4);
    cin >> N;
    A = new float* [N];
    for (int i = 0; i < N; i++) {
        A[i] = new float[N];
    }
    double time = 0.0;


    QueryPerformanceFrequency(&freq_);
    init();
    QueryPerformanceCounter(&begin_time);
    for (int k = 0; k < N; k++) {
        for (int j = k; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++) {
            for (int j = k + 1; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0.0;
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart - begin_time.QuadPart) / (double)(freq_.QuadPart);
    time = time * 1000;
    cout<<"chuan xing:  ";
    cout << time << endl;


    //直接调用opm
    init();
    QueryPerformanceCounter(&begin_time);
{
#pragma omp parallel for num_threads(THREADS_NUM)
        for (int k = 0; k < N; k++) {
            for (int j = k; j < N; j++) {
                A[k][j] = A[k][j] / A[k][k];
            }
            A[k][k] = 1.0;
            for (int i = k + 1; i < N; i++) {
                for (int j = k + 1; j < N; j++) {
                    A[i][j] = A[i][j] - A[i][k] * A[k][j];
                }
                A[i][k] = 0.0;
            }
        }
}
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart - begin_time.QuadPart) / (double)(freq_.QuadPart);
    time = time * 1000;
    cout<<"opm1 :  ";
    cout << time << endl;



    //均分
    int m = 0;
    init();
    QueryPerformanceCounter(&begin_time);
#pragma omp parallel num_threads(THREADS_NUM) shared(m)
    {
#pragma omp barrier
        int r = omp_get_thread_num();
        for (int k = 0; k < N; k++) {

            for (int i = k + 1; i < N; i++) {

                if ((i % THREADS_NUM) == r) {

                    A[i][k] = A[i][k] / A[k][k];

                    for (int j = k + 1; j < N; j++) {
                        A[i][j] = A[i][j] - A[i][k] * A[k][j];
                    }

                }
            }
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart - begin_time.QuadPart) / (double)(freq_.QuadPart);
    time = time * 1000;
    cout<<"opm2  :  ";
    cout << time << endl;



    //静态
    init();
    QueryPerformanceCounter(&begin_time);

    {

            int r = omp_get_thread_num();
#pragma omp parallel for num_threads(THREADS_NUM) schedule(static,(N/THREADS_NUM))
            for (int k = 0; k < N; k++) {

                for (int i = k + 1; i < N; i++) {

                    if ((i % THREADS_NUM) == r) {

                        A[i][k] = A[i][k] / A[k][k];

                        for (int j = k + 1; j < N; j++) {
                            A[i][j] = A[i][j] - A[i][k] * A[k][j];
                        }

                    }
                }
            }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart - begin_time.QuadPart) / (double)(freq_.QuadPart);
    time = time * 1000;
    cout <<"静态:   "<<time << endl;


    //动态
    init();
    QueryPerformanceCounter(&begin_time);
    {

            int r = omp_get_thread_num();
#pragma omp parallel for num_threads(THREADS_NUM) schedule(dynamic,64)
            for (int k = 0; k < N; k++) {

                for (int i = k + 1; i < N; i++) {

                    if ((i % THREADS_NUM) == r) {

                        A[i][k] = A[i][k] / A[k][k];

                        for (int j = k + 1; j < N; j++) {
                            A[i][j] = A[i][j] - A[i][k] * A[k][j];
                        }

                    }
                }
            }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart - begin_time.QuadPart) / (double)(freq_.QuadPart);
    time = time * 1000;
    cout <<"动态64:   "<< time << endl;


    //动态
    m = 0;
    init();
    QueryPerformanceCounter(&begin_time);
    {

        int r = omp_get_thread_num();
#pragma omp parallel for num_threads(THREADS_NUM) schedule(guided,16)
        for (int k = 0; k < N; k++) {

            for (int i = k + 1; i < N; i++) {

                if ((i % THREADS_NUM) == r) {

                    A[i][k] = A[i][k] / A[k][k];

                    for (int j = k + 1; j < N; j++) {
                        A[i][j] = A[i][j] - A[i][k] * A[k][j];
                    }

                }
            }
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart - begin_time.QuadPart) / (double)(freq_.QuadPart);
    time = time * 1000;
    cout << "动态16:   " << time << endl;


    //128位，simd加速+负载均衡
#pragma omp parallel num_threads(THREADS_NUM) shared(m)
    {
#pragma omp barrier
        int r = omp_get_thread_num();
        __m128 t1_1, t2, t3;

        for (int k = 0; k < N; k++) {

            for (int i = k + 1; i < N; i++) {

                if ((i % THREADS_NUM) == r) {

                    A[i][k] = A[i][k] / A[k][k];

                    int offset = (N - k - 1) % 4;
                    for (int j = k + 1; j < k + 1 + offset; j++) {
                        A[i][j] = A[i][j] - A[i][k] * A[k][j];
                    }
                    t2 = _mm_set_ps(A[i][k], A[i][k], A[i][k], A[i][k]);
                    for (int j = k + 1 + offset; j < N; j += 4) {
                        t3 = _mm_load_ps(A[k] + j);
                        t1_1 = _mm_load_ps(A[i] + j);
                        t2 = _mm_mul_ps(t2, t3);
                        t1_1 = _mm_sub_ps(t1_1, t2);
                        _mm_store_ps(A[i] + j, t1_1);
                    }
                }
            }
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart - begin_time.QuadPart) / (double)(freq_.QuadPart);
    time = time * 1000;
    cout<<"simd 128:  ";
    cout << time << endl;

/*
    //256位，simd加速+均分
    init();
#pragma omp parallel num_threads(THREADS_NUM) shared(m)
    {
#pragma omp barrier
        int r = omp_get_thread_num();
        __m256 t1_1, t2, t3;

        for (int k = 0; k < N; k++) {

            for (int i = k + 1; i < N; i++) {

                if ((i % THREADS_NUM) == r) {

                    A[i][k] = A[i][k] / A[k][k];

                    int offset = (N - k - 1) % 8;
                    for (int j = k + 1; j < k + 1 + offset; j++) {
                        A[i][j] = A[i][j] - A[i][k] * A[k][j];
                    }
                    float temp1[8] = { A[i][k],A[i][k],A[i][k],A[i][k],A[i][k],A[i][k],A[i][k],A[i][k] };
                    t2 = _mm256_loadu_ps(temp1);
                    for (int j = k + 1 + offset; j < N; j += 8) {
                        t3 = _mm256_load_ps(A[i] + j);
                        t1_1 = _mm256_load_ps(A[i] + j);
                        t2 = _mm256_mul_ps(t2, t3);
                        t1_1 = _mm256_sub_ps(t1_1, t2);
                        _mm256_store_ps(A[i] + j, t1_1);
                    }

                }
            }
        }
    }
    QueryPerformanceCounter(&end_time);
    time = (double)(end_time.QuadPart - begin_time.QuadPart) / (double)(freq_.QuadPart);
    time = time * 1000;
    cout<<"simd 256:  ";
    cout << time << endl;
*/
    return 0;

}
