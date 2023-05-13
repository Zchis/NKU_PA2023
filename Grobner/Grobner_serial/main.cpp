#include<fstream>
#include<ctime>
#include<Windows.h>
#include <immintrin.h>
#include <intrin.h>
#include<cerrno>
#include<string>
#include<sstream>
#include <emmintrin.h>
#include <immintrin.h>
#include<iostream>
using namespace std;
/*
1:  562   170    53
2:  562   170   53
3:  562   170   53
4:  1011  539   263
5:  2362  11706  453
6:  3799  2759  1953
7:  53399  6375  4535
53:  23075 1537453 14325
9:  39060 23904 14921
10: 43577 39477 54274
11: 535401 5724  756
*/

double time1 = 0, time2 = 0, time3 = 0, time4 = 0, time5 = 0, time6 = 0, time7 = 0;

using namespace std;
void serial() {
	ifstream ifs("E:\\GrobnrtGE\\data\\set3\\elimkey.txt");//�ļ��Ķ���
	ifstream if1("E:\\GrobnrtGE\\data\\set3\\elimtar.txt");
	ofstream out("E:\\GrobnrtGE\\data\\set3\\result.txt");

	int forline = 53;//����Ԫ��53��
	int line = 170;//��Ԫ��170��
	int artix = 562;//ÿһ����562��Ԫ��

	int xiaoyuanzi[562][562];
	for (int i = 0; i < 562; i++) {
		for (int j = 0; j < 562; j++) {
			xiaoyuanzi[i][j] = 0;
		}
	}



	int xiaoyuanhang[53][562];
	for (int i = 0; i < 53; i++) {
		for (int j = 0; j < 562; j++) {
			xiaoyuanhang[i][j] = 0;
		}
	}
	//��Ԫ�������



	int Result[562][562];
	for (int i = 0; i < 562; i++) {
		for (int j = 0; j < 562; j++) {
			Result[i][j] = 0;
		}
	}

	string s;
	int temp1 = 0;
	for (int i = 1; i <= 170; i++) {
		getline(if1, s);
		stringstream ss(s);
		int num;
		while (ss >> num) {
			//��ʼ��Ԫ
			if (temp1 == 0) {
				temp1 = num;
			}
			xiaoyuanzi[temp1][num] = 1;//ÿһ�е�һ��Ԫ��Ϊ���б�ţ�ͬһ����Ӧλ��Ϊ1
		}
		temp1 = 0;
	}
	//��Ԫ�Ӷ���


	string s1;
	for (int i = 1; i <= 53; i++) {
		getline(ifs, s1);
		stringstream ss(s1);
		int num;
		while (ss >> num) {
			//��ʼ��Ԫ
			xiaoyuanhang[i - 1][num] = 1;
		}
	}
	//��Ԫ�ӵĶ���


	long long head3, tail3, freq3; // timers
	// similar to CLOCKS_PER_SEC
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq3);
	// start time
	QueryPerformanceCounter((LARGE_INTEGER*)&head3);

	int t = 0;//������Ԫ
	int T = 0;
	while (true) {
		if (t == 53) {
			break;
		}
		for (int j = 129; j >= 0; j--) {
			if (xiaoyuanhang[t][j] == 1) {//����1
				T = 1;
				if (xiaoyuanzi[j][j] == 0) {//������Ԫ��
					for (int k = 0; k < 562; k++) {
						xiaoyuanzi[j][k] = xiaoyuanhang[t][k];
						Result[j][k] = xiaoyuanhang[t][k];
					}
					t++;
					//cout << "easeasd" << result[j][j] << " " << j << endl;
					//cout << t << endl;
					break;
				}
				else {
					for (int k = 0; k < 562; k++) {
						if (xiaoyuanzi[j][k] == xiaoyuanhang[t][k]) {
							xiaoyuanhang[t][k] = 0;
						}
						else {
							xiaoyuanhang[t][k] = 1;
						}
					}
				}
			}
		}
		if (T == 1) {
			T = 0;
		}
		else {
			//t=0��ȫ0������Ԫ����
			t++;
		}
	}

	QueryPerformanceCounter((LARGE_INTEGER*)&tail3);
	time1 += (tail3 - head3) * 1000.0 / freq3;
	cout << "time:    " << (tail3 - head3) * 1000.0 / freq3 << "ms" << endl;
	/*
	int temp2 = 0;
	for (int i = 129; i >= 0; i--)
	{
		for (int j = 129; j >= 0; j--) {
			if (Result[i][j] == 1) {
				cout << j << " ";
				temp2 = 1;
			}
			if (j == 0 && temp2 == 1) {
				cout << endl;
				temp2 = 0;
			}
		}
	}*/

	ifs.close();
}
int main()
{
    serial();
    return 0;
}
