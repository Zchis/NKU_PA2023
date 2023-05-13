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

#include<immintrin.h>  // AVX

#include<pthread.h>
#include<semaphore.h>
#include<mutex>
#include<condition_variable>

using namespace std;

#define THREAD_NUM 8  // �߳�����


//--------------------------------------��������--------------------------------------


int n = 0;  // �����С
const int k = 1;

int ek_num;  // ������Ԫ�Ӹ���
int et_num;  // ���뱻��Ԫ������

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

�ļ��ڲ�������

MyBitSet���ݽṹ����Ӧ�ô����ݽṹ�Ĵ����㷨�������ţ��Լ�����PThread�㷨��SIMD�㷨���Լ��ۺ�
*/

//------------------------------------���ݵ��빤��------------------------------------
const int column_num_c = 8399;
const int ek_num_c = 6375;  // ������Ԫ�Ӹ���
const int et_num_c = 4535;  // ���뱻��Ԫ������

string dir = "E:\\GrobnrtGE\\data\\set7\\";
stringstream ss;


int bit_size = column_num_c / 32 + 1;
class MyBitSet {
public:
	int head;  // ����
	int* content;

	MyBitSet() {
		head = -1; content = new int[bit_size];
		for (int i = 0; i < bit_size; i++) content[i] = 0;
	}

	// bool operator[](size_t index) {}

	MyBitSet& operator^=(const MyBitSet& b) {  // Ĭ����������bitset������ͬ
		for (int i = 0; i < bit_size; i++) content[i] ^= b.content[i];
		for (int i = 0; i < bit_size; i++)
		{
			for (int j = 0; j < 32; j++)
			{
				if ((content[i] & (1 << j)))
				{
					head = i * 32 + j;
					return *this;
				}
			}
		}
		head = -1;
		return *this;
	}

	MyBitSet& my_xor_AVX(const MyBitSet& b) {
		__m256i v_this, v_b;
		int i = 0;
		for (i; i < bit_size - 8; i += 8) {
			v_this = _mm256_loadu_si256((__m256i*) & content[i]);
			v_b = _mm256_loadu_si256((__m256i*) & b.content[i]);
			v_this = _mm256_xor_si256(v_this, v_b);
			_mm256_storeu_si256((__m256i*) & content[i], v_this);
		}
		for (i; i < bit_size; i++)
		{
			content[i] ^= b.content[i];
		}
		for (int i = 0; i < bit_size; i++)
		{
			for (int j = 0; j < 32; j++)
			{
				if ((content[i] & (1 << j)))
				{
					head = i * 32 + j;
					return *this;
				}
			}
		}
		head = -1;
		return *this;
	}


	int test(int index) {
		return content[index / 32] & (1 << (index % 32)) ? 1 : 0;  // Ѱַ��ʽ
	}

	void set(int index) {  // ��λ
		content[index / 32] |= (1 << (index % 32));
	}

	bool any() {
		for (int i = 0; i < bit_size; i++) if (content[i]) return true;
		return false;
	}

private:

};


bitset<column_num_c> eks_c[column_num_c];  // ��Ԫ�ӣ�����һЩ���ڼ���������
bitset<column_num_c> ets_c[et_num_c];

int lp_ets_c[et_num_c];
int lp_eks_c[column_num_c];


MyBitSet eks[column_num_c];
MyBitSet ets[et_num_c];



long long head, tail, freq;

//------------------------------------�߳����ݽṹ------------------------------------

typedef struct {
	int t_id;  // �߳� id
	int tasknum;  // ��������
}PT_EliminationParam;

typedef struct {
	int t_id; //�߳� id
	size_t index;
}PT_GetHeadParam;

//-------------------------------------�ź�������-------------------------------------
sem_t sem_gethead;
sem_t sem_workerstart[THREAD_NUM]; // ÿ���߳����Լ�ר�����ź���
sem_t sem_workerend[THREAD_NUM];

sem_t sem_leader;
sem_t sem_Divsion[THREAD_NUM - 1];
sem_t sem_Elimination[THREAD_NUM - 1];

//------------------------------------barrier����-------------------------------------
pthread_barrier_t barrier_gethead;
pthread_barrier_t _Barrier_Elimination;
pthread_barrier_t _Barrier_0_Operation;

//-------------------------------���������ͻ���������---------------------------------
pthread_mutex_t _mutex;
pthread_cond_t _cond;

//--------------------------------------�̺߳���--------------------------------------


unordered_set<int> ek_set;
void* PT_Elimination_Static(void* param) {
	PT_EliminationParam* tt = (PT_EliminationParam*)param;
	int t_id = tt->t_id;
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{
			pthread_barrier_wait(&_Barrier_Elimination);  // ��֤ͬ��
			if (t_id == 0)  // ����̸߳���������Ԫ�Ӳ���
			{
				for (size_t j = 0; j < et_num_c; j++)
				{
					if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
					{
						eks[i] = ets[j];
						ets[j].head = -1;
						break;
					}
				}
			}
			pthread_barrier_wait(&_Barrier_0_Operation);  // ��0���߳̽�������Ԫ���������֮����Ҫ��֤����ͬ��
		}
		for (size_t j = t_id; j < et_num_c; j += THREAD_NUM)  // ѭ�����ֲ��л�
		{
			if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets[j] ^= eks[i];
				for (int k = i; k < column_num_c; k++)
				{
					if (ets[j].test(k))
					{
						ets[j].head = k;
						break;
					}
					if (k == column_num_c - 1) ets[j].head = -1;
				}
			}
		}
	}
	pthread_exit(nullptr);
	return nullptr;
}

void* PT_Elimination_Static_Block(void* param) {
	PT_EliminationParam* tt = (PT_EliminationParam*)param;
	int t_id = tt->t_id;
	int tasknum = t_id == THREAD_NUM - 1 ? et_num_c % (THREAD_NUM - 1) : et_num_c / (THREAD_NUM - 1);
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{
			pthread_barrier_wait(&_Barrier_Elimination);  // ��֤ͬ��
			if (t_id == 0)  // ����̸߳���������Ԫ�Ӳ���
			{
				for (size_t j = 0; j < et_num_c; j++)
				{
					if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
					{
						eks[i] = ets[j];
						ets[j].head = -1;
						break;
					}
				}
			}
			pthread_barrier_wait(&_Barrier_0_Operation);  // ��0���߳̽�������Ԫ���������֮����Ҫ��֤����ͬ��
		}
		for (size_t j = t_id * tasknum, k = 0; k < tasknum; j++, k++)  // �黮�ֲ��л�
		{
			if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets[j] ^= eks[i];
				for (int k = i; k < column_num_c; k++)
				{
					if (ets[j].test(k))
					{
						ets[j].head = k;
						break;
					}
					if (k == column_num_c - 1) ets[j].head = -1;
				}
			}
		}
	}
	pthread_exit(nullptr);
	return nullptr;
}


void* PT_Elimination_Static_AVX(void* param) {
	PT_EliminationParam* tt = (PT_EliminationParam*)param;
	int t_id = tt->t_id;
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{
			pthread_barrier_wait(&_Barrier_Elimination);  // ��֤ͬ��
			if (t_id == 0)  // ����̸߳���������Ԫ�Ӳ���
			{
				for (size_t j = 0; j < et_num_c; j++)
				{
					if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
					{
						eks[i] = ets[j];
						ets[j].head = -1;
						break;
					}
				}
			}
			pthread_barrier_wait(&_Barrier_0_Operation);  // ��0���߳̽�������Ԫ���������֮����Ҫ��֤����ͬ��
		}
		for (size_t j = t_id; j < et_num_c; j += THREAD_NUM)  // ѭ�����ֲ��л�
		{
			if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets[j].my_xor_AVX(eks[i]);
				for (int k = i; k < column_num_c; k++)
				{
					if (ets[j].test(k))
					{
						ets[j].head = k;
						break;
					}
					if (k == column_num_c - 1) ets[j].head = -1;
				}
			}
		}
	}
	pthread_exit(nullptr);
	return nullptr;
}

void* PT_Elimination_Static_Block_AVX(void* param) {
	PT_EliminationParam* tt = (PT_EliminationParam*)param;
	int t_id = tt->t_id;
	int tasknum = t_id == THREAD_NUM - 1 ? et_num_c % (THREAD_NUM - 1) : et_num_c / (THREAD_NUM - 1);
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{
			pthread_barrier_wait(&_Barrier_Elimination);  // ��֤ͬ��
			if (t_id == 0)  // ����̸߳���������Ԫ�Ӳ���
			{
				for (size_t j = 0; j < et_num_c; j++)
				{
					if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
					{
						eks[i] = ets[j];
						ets[j].head = -1;
						break;
					}
				}
			}
			pthread_barrier_wait(&_Barrier_0_Operation);  // ��0���߳̽�������Ԫ���������֮����Ҫ��֤����ͬ��
		}
		for (size_t j = t_id * tasknum, k = 0; k < tasknum; j++, k++)  // �黮�ֲ��л�
		{
			if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets[j].my_xor_AVX(eks[i]);
				for (int k = i; k < column_num_c; k++)
				{
					if (ets[j].test(k))
					{
						ets[j].head = k;
						break;
					}
					if (k == column_num_c - 1) ets[j].head = -1;
				}
			}
		}
	}
	pthread_exit(nullptr);
	return nullptr;
}


//-------------------------------------PThread�㷨------------------------------------

void PT_GrobnerGE_Static()
{
	pthread_barrier_init(&_Barrier_Elimination, nullptr, THREAD_NUM);
	pthread_barrier_init(&_Barrier_0_Operation, nullptr, THREAD_NUM);
	//int tasknum = eks.size() / THREAD_NUM;
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_EliminationParam* param = (PT_EliminationParam*)malloc(THREAD_NUM * sizeof(PT_EliminationParam));  // �洢�̲߳���
	for (int t_id = 0; t_id < THREAD_NUM; t_id++) {
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], nullptr, PT_Elimination_Static, (void*)&param[t_id]);
	}

	for (int t_id = 0; t_id < THREAD_NUM; t_id++) pthread_join(handles[t_id], nullptr);

	pthread_barrier_destroy(&_Barrier_Elimination);
	pthread_barrier_destroy(&_Barrier_0_Operation);
}

void PT_GrobnerGE_Static_Block()
{
	pthread_barrier_init(&_Barrier_Elimination, nullptr, THREAD_NUM);
	pthread_barrier_init(&_Barrier_0_Operation, nullptr, THREAD_NUM);
	//int tasknum = eks.size() / THREAD_NUM;
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_EliminationParam* param = (PT_EliminationParam*)malloc(THREAD_NUM * sizeof(PT_EliminationParam));  // �洢�̲߳���
	for (int t_id = 0; t_id < THREAD_NUM; t_id++) {
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], nullptr, PT_Elimination_Static_Block, (void*)&param[t_id]);
	}

	for (int t_id = 0; t_id < THREAD_NUM; t_id++) pthread_join(handles[t_id], nullptr);

	pthread_barrier_destroy(&_Barrier_Elimination);
	pthread_barrier_destroy(&_Barrier_0_Operation);
}

void PT_GrobnerGE_Static_AVX()
{
	pthread_barrier_init(&_Barrier_Elimination, nullptr, THREAD_NUM);
	pthread_barrier_init(&_Barrier_0_Operation, nullptr, THREAD_NUM);
	//int tasknum = eks.size() / THREAD_NUM;
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_EliminationParam* param = (PT_EliminationParam*)malloc(THREAD_NUM * sizeof(PT_EliminationParam));  // �洢�̲߳���
	for (int t_id = 0; t_id < THREAD_NUM; t_id++) {
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], nullptr, PT_Elimination_Static_AVX, (void*)&param[t_id]);
	}

	for (int t_id = 0; t_id < THREAD_NUM; t_id++) pthread_join(handles[t_id], nullptr);

	pthread_barrier_destroy(&_Barrier_Elimination);
	pthread_barrier_destroy(&_Barrier_0_Operation);
}

void PT_GrobnerGE_Static_Block_AVX()
{
	pthread_barrier_init(&_Barrier_Elimination, nullptr, THREAD_NUM);
	pthread_barrier_init(&_Barrier_0_Operation, nullptr, THREAD_NUM);
	//int tasknum = eks.size() / THREAD_NUM;
	pthread_t* handles = (pthread_t*)malloc(THREAD_NUM * sizeof(pthread_t));  // Ϊ�߳̾�������ڴ�ռ�
	PT_EliminationParam* param = (PT_EliminationParam*)malloc(THREAD_NUM * sizeof(PT_EliminationParam));  // �洢�̲߳���
	for (int t_id = 0; t_id < THREAD_NUM; t_id++) {
		param[t_id].t_id = t_id;
		pthread_create(&handles[t_id], nullptr, PT_Elimination_Static_Block_AVX, (void*)&param[t_id]);
	}

	for (int t_id = 0; t_id < THREAD_NUM; t_id++) pthread_join(handles[t_id], nullptr);

	pthread_barrier_destroy(&_Barrier_Elimination);
	pthread_barrier_destroy(&_Barrier_0_Operation);
}



//------------------------------------���㸨������------------------------------------

int find_first_bitset(const bitset<column_num_c>& b) {
	for (int i = 0; i < column_num_c; i++) if (b.test(i)) return i;
	return -1;
}

int get_head_lp_c(size_t index) {  // ��ȡ����
	for (int i = 0; i < column_num_c; i++)
	{
		cout << index << " " << i << endl;
		if (index == lp_eks_c[i]) return i;
		else if (lp_eks_c[i] == -1) break;
	}
	return -1;
}

void xor_static_AVX(bitset<column_num_c>& b1, const bitset<column_num_c>& b2) {  // b1 ^= b2����֤b1��b2�ȳ�
	for (size_t _Wpos = 0; _Wpos < column_num_c; ++_Wpos) {
		//b1.set(_Wpos, b1._Subscript(_Wpos) ^ b2._Subscript(_Wpos));
	}
	//b1 ^= b2;
}

//--------------------------------------�㷨����--------------------------------------

void GrobnerGE()
{
	int lp = -1;
	for (int i = 0; i < et_num_c; i++)
	{
		while (ets[i].any())
		{
			lp = ets[i].head;
			if (eks[lp].test(lp))  // ˵�����ڶ�Ӧ��Ԫ��
			{
				ets[i] ^= eks[lp];
			}
			else
			{
				eks[lp] = ets[i];
				ets[i].head = -1;
				break;
			}
		}
	}
}

void GrobnerGE_LP_static()  // ��̬λ���µĸ�˹��ȥ
{
	int lp = -1;
	for (int i = 0; i < et_num_c; i++)
	{
		while (ets_c[i].any())
		{
			lp = lp_ets_c[i];
			if (eks_c[lp].test(lp))  // ˵�����ڶ�Ӧ��Ԫ��
			{
				ets_c[i] ^= eks_c[lp];
				lp_ets_c[i] = find_first_bitset(ets_c[i]);
			}
			else
			{
				eks_c[lp] = ets_c[i];
				lp_ets_c[i] = -1;
				break;
			}
		}
	}
}

void GrobnerGE_LP_static_headopt()  // ��̬λ���µĸ�˹��ȥ��ѡȡ�����һ���Ż�
{
	int lp = -1;
	for (int i = 0; i < et_num_c; i++)
	{
		while (ets_c[i].any())
		{
			lp = lp_ets_c[i];
			if (eks_c[lp].test(lp))  // ˵�����ڶ�Ӧ��Ԫ��
			{
				ets_c[i] ^= eks_c[lp];
				for (int j = lp; j < column_num_c; j++)
					if (ets_c[i].test(j)) {
						lp_ets_c[i] = j;
						break;
					}
			}
			else
			{
				eks_c[lp] = ets_c[i];
				lp_ets_c[i] = -1;
				break;
			}
		}
	}
}

void GrobnerGE_Rearrange() {
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{  // ����̸߳���������Ԫ�Ӳ���
			for (size_t j = 0; j < et_num_c; j++)
			{
				if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
				{
					eks[i] = ets[j];
					ets[j].head = -1;
					break;
				}
			}
		}
		for (size_t j = 0; j < et_num_c; j++)  // ѭ�����ֲ��л�
		{
			if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets[j] ^= eks[i];
				for (int k = i; k < column_num_c; k++)
				{
					if (ets[j].test(k))
					{
						ets[j].head = k;
						break;
					}
					if (k == column_num_c - 1) ets[j].head = -1;
				}
			}
		}
	}
}

void GrobnerGE_Rearrange_AVX() {
	for (int i = 0; i < column_num_c; i++)  // ȡÿ����Ԫ�ӣ��Ա���Ԫ�н��в��������ڲ��л�
	{
		if (!eks[i].test(i))  // ��Ԫ�ӱ������ʼ��ʱ���㡰�кš� = �����������
		{  // ����̸߳���������Ԫ�Ӳ���
			for (size_t j = 0; j < et_num_c; j++)
			{
				if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
				{
					eks[i] = ets[j];
					ets[j].head = -1;
					break;
				}
			}
		}
		for (size_t j = 0; j < et_num_c; j++)  // ѭ�����ֲ��л�
		{
			if (i == ets[j].head)  // ˵�����ڶ�Ӧ����Ԫ��
			{
				ets[j].my_xor_AVX(eks[i]);
				//ets[j] ^= eks[i];
				for (int k = i; k < column_num_c; k++)
				{
					if (ets[j].test(k))
					{
						ets[j].head = k;
						break;
					}
					if (k == column_num_c - 1) ets[j].head = -1;
				}
			}
		}
	}
}

//------------------------------------������Ժ���------------------------------------

void reverse_output_c()
{
	ofstream outp(dir + "output_c.txt");
	for (int i = 0; i < et_num_c; i++)
	{
		for (int j = 0; j < n; j++) if (ets_c[i].test(j)) outp << n - j - 1 << " ";
		outp << endl;
	}
	outp.close();
}

void reverse_output_MyB()
{
	ofstream outp(dir + "output_MyB.txt");
	for (int i = 0; i < et_num_c; i++)
	{
		for (int j = 0; j < n; j++) if (ets[i].test(j)) outp << n - j - 1 << " ";
		outp << endl;
	}
	outp.close();
}

//------------------------------------���ݶ�ȡ����------------------------------------

void readData_reverse_bitset()
{  // �����������
	string inek, inet;
	stringstream ss_inek, ss_inet;
	ifstream inElimKey(dir + "elimkey.txt");  // ��Ԫ��
	ifstream inElimTar(dir + "elimtar.txt");  // ����Ԫ��
	int inek_loc, p_ek = 0, inet_loc, p_et = 0;  // �������ݶ���


	while (true)  // ��ȡ��Ԫ��
	{
		getline(inElimKey, inek);
		ss_inek = stringstream(inek);
		while (ss_inek >> inek)
		{
			inek_loc = stoi(inek);
			cout << inek_loc << " ";
			eks[p_ek].set(n - inek_loc - 1);
		}
		if (inek.empty())
		{
			cout << "ek_complete" << endl;
			break;
		}
		//cout << eks[p_ek] << endl;
		p_ek++;
	}
	while (true)  // ��ȡ����Ԫ��
	{
		getline(inElimTar, inet);
		ss_inet = stringstream(inet);
		while (ss_inet >> inet)
		{
			inet_loc = stoi(inet);
			//cout << inet_loc << " ";
			ets[p_et].set(n - inet_loc - 1);
		}
		if (inet.empty())
		{
			//cout << "et_complete" << endl;
			break;
		}
		//cout << ets[p_et] << endl;
		p_et++;
	}
}

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
			if (lp == -1) lp = n - inek_loc - 1;
			//cout << inek_loc << " ";
			eks_c[lp].set(n - inek_loc - 1);
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
				lp = n - inet_loc - 1;
				lp_ets_c[p_et] = lp;
			}
			//cout << inet_loc << " ";
			ets_c[p_et].set(n - inet_loc - 1);
		}
		if (inet.empty()) break;
		lp = -1;
		p_et++;
		//cout << ets_c[p_et] << endl;
	}
	//cout << "et_complete" << endl;
	inElimKey.close();
	inElimTar.close();
	cout << "init_complete" << endl;
}

void readData_reverse_MyB()
{  // ����������ݣ����뾲̬λ��
    cout<<"begin!"<<endl;
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
			if (lp == -1)
			{
				lp = column_num_c - inek_loc - 1;
				eks[lp].head = lp;
			}
			//cout << inek_loc << " ";
			eks[lp].set(column_num_c - inek_loc - 1);

		}
		lp = -1;  p_ek++;
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
				ets[p_et].head = lp;
			}
			//cout << inet_loc << " ";
			ets[p_et].set(column_num_c - inet_loc - 1);
		}
		lp = -1;  p_et++;
		if (inet.empty()) break;
		//cout << ets_c[p_et] << endl;
	}
	//cout << "et_complete" << endl;
	inElimKey.close();
	inElimTar.close();
	//cout << "init_complete" << endl;

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

void init_MyB() {
	for (int i = 0; i < column_num_c; i++)
	{

	}
	readData_reverse_MyB();  // �����ʼ����Ԫ�Ӻͱ���Ԫ������
	cout << "init_complete" << endl;
}

int main() {
	ifstream inParam(dir + "param.txt");
	inParam >> n;  // ��������С
	inParam >> ek_num;  // ���������Ԫ�Ӹ���
	inParam >> et_num;  // ���뱻��Ԫ������
	//��̬�㷨����Ҫ�ֶ����ĳ����е�ȫ�ֳ���

	cout << "  n   " << n << "   ek_num geshu   " << ek_num << "   et_num hang " << et_num << endl;

	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	//-----------------------------------------------------------------

	init_MyB();


	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE();  // ��̬λ���洢�ľ���������˹��ȥ
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE: " << (tail - head) * 1000.0 / freq
		<< " ms" << endl;
	//reverse_output_MyB();

	//-----------------------------------------------------------------

	//init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_Rearrange();  // �����ݽṹ�������㷨
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_Rearrange: " << (tail - head) * 1000.0 / freq
		<< " ms" << endl;
	//reverse_output_MyB();

	//-----------------------------------------------------------------

	//init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	GrobnerGE_Rearrange_AVX();  // �����ݽṹ�������㷨+SIMD
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "GrobnerGE_Rearrange_AVX: " << (tail - head) * 1000.0 / freq
		<< " ms" << endl;
	//reverse_output_MyB();


	//-----------------------------------------------------------------

	//init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static();  // GGELSH���߳�ѭ�����ְ汾
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static: " << (tail - head) * 1000.0 / freq
		<< " ms" << endl;
	//reverse_output_MyB();

	//-----------------------------------------------------------------

	//init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static_Block();  // GGELSH���߳̿黮�ְ汾
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static_Block: " << (tail - head) * 1000.0 / freq
		<< " ms" << endl;
	//reverse_output_MyB();

	//-----------------------------------------------------------------

	//init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static_AVX();  // GGELSH���߳̿黮�ְ汾+SIMD
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static_AVX: " << (tail - head) * 1000.0 / freq
		<< " ms" << endl;
	reverse_output_MyB();

	//-----------------------------------------------------------------

	//init_MyB();
	QueryPerformanceCounter((LARGE_INTEGER*)&head);
	PT_GrobnerGE_Static_Block_AVX();  // GGELSH���߳̿黮�ְ汾+SIMD
	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
	cout << "PT_GrobnerGE_Static_Block_AVX: " << (tail - head) * 1000.0 / freq
		<< " ms" << endl;
	//reverse_output_MyB();



}


