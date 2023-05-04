#include<iostream>
#include<immintrin.h>
#include<random>
#include<stdio.h>
#include<stdio.h>
#include<time.h>
using namespace std;
int main() {
	const int n = 2000;
	float A[n][n];
	for (int i = 0; i < n; i++) {
		A[i][i] = 1.0;
		for (int j = 0; j < n; j++)
		{
			A[i][j] = rand();
		}
	}
	for (int k = 0; k < n; k++) {
		for (int i = k + 1; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] += A[k][j];
			}
		}
	}

	struct timespec sts, ets;
	timespec_get(&sts, TIME_UTC);


	for (int k = 0; k < n; k++) {
		__m128 vt = _mm_set1_ps(A[k][k]);
		//float32x4_t vt = vmovq_n_f32(A[k][k]);//传值初始化
		for (int j = k + 1; j < n; j += 4) {//一定要是4的倍数
			//float32x4_t va = vld1q_f32(&A[k][j]);//传递地址
			__m128 va = _mm_loadu_ps(&A[k][j]);//不要求对齐的
			//__m128 va = _mm_load_ps(&A[k][j]);//要求16字节对齐
			//va = vdivq_f32(va, vt);
			va = _mm_div_ps(va, vt);
			//vst1q_f32(&A[k][j], va);
			_mm_store_ps(&A[k][j], va);
		}
		//省略一部
		for (int i = k + 1; i < n; i++) {
			//float32x4_t vaik = vmovq_n_f32(A[i][k]);
			__m128 vaik = _mm_set1_ps(A[i][k]);
			int j = k + 1;
			for (; j+4 <= n; j += 4) {
				//float32x4_t vakj = vld1q_f32(&A[k][j]);//加载
				__m128 vakj = _mm_loadu_ps(&A[k][j]);//不要求对齐
				//__m128 vakj = _mm_load_ps(&A[k][j]);//要求16字节对齐

				//float32x4_t vaij = vld1q_f32(&A[i][j]);
				__m128 vaij = _mm_loadu_ps(&A[i][j]);
				//__m128 vaij = _mm_load_ps(&A[i][j]);//要求16字节对齐

				//float32x4_t vx = vmovq_n_f32(0);//是否需要初始化
				__m128 vx = _mm_set1_ps(0);

				//vx = float32x4_t vmulq_f32(vakj, vaik);//乘法
				vx = _mm_mul_ps(vakj, vaik);//乘法
				vaij = _mm_sub_ps(vaij, vx);//减法
				_mm_store_ps(&A[i][j], vaij);//加载回内存
			}
			for (; j < n; j++) {
				A[i][j] = A[i][j] - A[k][j] * A[i][k];
			}
		}
	}
	timespec_get(&ets, TIME_UTC);
	time_t dsec = ets.tv_sec-sts.tv_sec;
	long dnsec = ets.tv_nsec-sts.tv_nsec;
	if (dnsec < 0) {
		dsec--;
		dnsec += 1000000000ll;
	}
	printf("%lld.%09llds\n", dsec, dnsec);
}