#include<iostream>
#include<random>
#include<bits/stdc++.h>
#include<sys/time.h>

using namespace std;
const int N = 1200;
const double eps = 1e-6;
int n;
double a[N][N];
struct timespec sts, ets;
int main()
{
	for (int i = 0; i < n; i++) {
		a[i][i] = 1.0;
		for (int j = 0; j < n; j++) {
			a[i][j] = rand();
		}
	}
	for (int k = 0; k < n; k++) {
		for (int i = k + 1; i < n; i++) {
			for (int j = 0; j < n; j++) {
				a[i][j] += a[k][j];
			}
		}
	}
	timespec_get(&sts, TIME_UTC);
	int r, c; 
	for (r = 0, c = 0; c < n; c++)
	{
		int t = r;
		for (int i = r; i < n; i++)
			if (fabs(a[i][c]) > fabs(a[t][c]))
				t = i;
		if (fabs(a[t][c]) < eps)continue;
		for (int i = c; i < n + 1; i++)swap(a[t][i], a[r][i]);
		for (int i = n; i >= c; i--)a[r][i] /= a[r][c];
		for (int i = r + 1; i < n; i++)
			if (fabs(a[i][c]) > eps)
				for (int j = n; j >= c; j--)
					a[i][j] -= a[r][j] * a[i][c];
		r++;
	}
	timespec_get(&ets, TIME_UTC);
	time_t dsec = ets.tv_sec - sts.tv_sec;
	long dnsec = ets.tv_nsec - sts.tv_nsec;

	if (dnsec < 0) {
		dsec--;
		dnsec += 1000000000ll;
	}
	printf("%ld.%09lds\n", dsec, dnsec);
	cout << "Over!" << endl;

}

