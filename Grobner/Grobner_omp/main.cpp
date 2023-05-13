#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<ctime>
#include<bitset>
#include<vector>
#include<string>
#include<sstream>

#include<unordered_set>

#include<bitset>

#include<windows.h>
#include<omp.h>


using namespace std;

#define THREAD_NUM 8  // �߳�����


int n = 0;  // �����С
const int k = 1;

/* �������ݹ�ģ��
1:  130   22    8
2:  254   106   53
3:  562   170   53
4:  1011  539   263
5:  2362  1226  453
6:  3799  2759  1953
7:  8399  6375  4535
8:  23075 18748 14325
9:  39060 23904 14921
10: 43577 39477 54274
11: 85401 5724  756


�ļ��ڲ�������std::bitset�汾�������Ŵ����omp�㷨
*/

const int column_num_c = 3799;
const int ek_num_c = 2759;  // ������Ԫ�Ӹ���
const int et_num_c = 1953;  // ���뱻��Ԫ������


bitset<column_num_c> eks_c[column_num_c];  // ��Ԫ�ӣ�����һЩ���ڼ���������
bitset<column_num_c> ets_c[et_num_c];

int lp_ets_c[et_num_c];
int lp_eks_c[column_num_c];

long long head, tail, freq;

//------------------------------------���ݵ��빤��------------------------------------

string dir = "E:\\GrobnrtGE\\data\\set6\\";
stringstream ss;

//------------------------------------���㸨������------------------------------------
void GrobnerGE_OMP()
{
    int i, j, k;
    #pragma omp parallel num_threads(THREAD_NUM), private(i, j, k)
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks_c[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{
            #pragma omp barrier
            #pragma omp single
            for (size_t j = 0; j < et_num_c; j++)
            {
                if (i == lp_ets_c[j])  // ˵�����ڶ�Ӧ����Ԫ��
                {
                    eks_c[i] = ets_c[j];
                    lp_ets_c[j] = -1;
                    break;
                }
            }
            //#pragma omp barrier
		}
        #pragma omp for schedule(guided)
		for (int j = 0; j < et_num_c; j ++)  // ѭ�����ֲ��л�
		{
			if (i == lp_ets_c[j])  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets_c[j] ^= eks_c[i];
				for (int k = i; k < column_num_c; k++)
				{
					if (ets_c[j].test(k))
					{
						lp_ets_c[j] = k;
						break;
					}
					if (k == column_num_c - 1) lp_ets_c[j] = -1;
				}
			}
		}
	}
}

int find_first_bitset(const bitset<column_num_c>& b) {
	for (int i = 0; i < column_num_c; i++) if (b.test(i)) return i;
	return -1;
}


//------------------------------------������Ժ���------------------------------------

void reverse_output_c()
{
	ofstream outp(dir + "output_c_omp.txt");
	for (int i = 0; i < et_num_c; i++)
	{
		for (int j = 0; j < column_num_c; j++) if (ets_c[i].test(j)) outp << column_num_c - j - 1 << " ";
		outp << endl;
	}
    cout<<"over   "<<endl;
	outp.close();
}

//------------------------------------���ݶ�ȡ����------------------------------------

void readData_reverse_bitset_c()
{  // ����������ݣ����뾲̬λ��
	string inek, inet;
	stringstream ss_inek, ss_inet;
	ifstream inElimKey(dir + "elimkey.txt");  // ��Ԫ��
	ifstream inElimTar(dir + "elimtar.txt");  // ����Ԫ��
	int inek_loc, p_ek = 0, inet_loc, p_et = 0;  // �������ݶ���
	int lp = -1;
	while (true)  // ��ȡ��Ԫ��
	{
		getline(inElimKey, inek);
		ss_inek = stringstream(inek);
		while (ss_inek >> inek)
		{
			inek_loc = stoi(inek);
			if (lp == -1) lp = column_num_c - inek_loc - 1;
			//cout << inek_loc << " ";
			eks_c[lp].set(column_num_c - inek_loc - 1);
		}
		lp = -1, p_ek++;
		if (inek.empty()) break;
		//cout << eks_c[p_ek] << endl;
	}
	//cout << "ek_complete" << endl;

	while (true)  // ��ȡ����Ԫ��
	{
		getline(inElimTar, inet);
		ss_inet = stringstream(inet);
		while (ss_inet >> inet)
		{
			inet_loc = stoi(inet);
			if (lp == -1)
			{
				lp = column_num_c - inet_loc - 1;
				lp_ets_c[p_et] = lp;
			}
			//cout << inet_loc << " ";
			ets_c[p_et].set(column_num_c - inet_loc - 1);
		}
		if (inet.empty()) break;
		lp = -1;
		p_et++;
		//cout << ets_c[p_et] << endl;
	}
	//cout << "et_complete" << endl;
	inElimKey.close();
	inElimTar.close();
}

void init_c() {
	for (int i = 0; i < column_num_c; i++)  // ��ʼ����������һ���㷨�еĲ�������
	{
		eks_c[i] = *(new bitset<column_num_c>);
		lp_eks_c[i] = -1;
	}
	readData_reverse_bitset_c();  // �����ʼ����Ԫ�Ӻͱ���Ԫ������
	cout << "init_complete" << endl;
}





int main() {
	//��̬�㷨����Ҫ�ֶ����ĳ����е�ȫ�ֳ���

	cout << "  n   " << column_num_c << ", ek_num_c   " << ek_num_c << ", et_num_c  " << et_num_c << endl;

	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	//-----------------------------------------------------------------

	init_c();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_OMP();  // ��̬λ���洢�ľ���������˹��ȥ
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_OMP: " << (tail - head) * 1000.0 / freq
		<< " ms" << endl;
    reverse_output_c();

	//-----------------------------------------------------------------
}


