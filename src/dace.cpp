/*
 *   DACE: A Scalable DP-means Algorithm for Clustering Extremely Large Sequence Data 
 *   Author: Linghao Jiang, Yichao Dong, Ning Chen, Ting Chen
 *   
 *   URL: https://github.com/tinglab/DACE
 * 
 *   program written by 
 *        Linhao Jiang, jianglh.thu@gmail.com, Tsinghua University
 *                   at
 *        Ting Chen's lab, tingchen@mail.tsinghua.edu.cn, Tsinghua University
 * 
 */

#include "mpi.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <unistd.h>
#include <cstdio>
#include <cmath>
#include <list>
#include <map>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <vector>
#include <pthread.h>
#include <limits.h>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "cls_common.h"
#include "cls_log.h"
#include "cls_alignment.h"
#include "cls_suffix_array.h"
#include "cls_cluster.h"
#include "cls_kmer.h"

#define MAX_SEQ_LEN 10000
#define MPI_TAG_SENDDATA 2001000
#define MPI_TAG_SENDDATA_LEN 2002000
#define MPI_TAG_RESULT 2003000
#define MPI_TAG_RESULT_LEN 2004000

#define PCA_KMER 4
#define PCA_U_LEN (1 << (PCA_KMER * 2))
#define PCA_SELECT 500

#define MPI_TAG_REQ_BLOCKID 0
#define MPI_TAG_ANS_BLOCKID 0

#define FIX_RAND_SEED

#ifdef USING_PCA
#include "f2c.h"
#include "clapack.h"
#endif

using namespace std;

FILE *fp;

char	*INPUT_FILE,	 /* input file name */
	*OUTPUT_FILE;	 /* output file name */

int	MPI_RANK,	 /* MPI rank */
	MPI_SIZE;	 /* MPI process number */

long	LEVEL = 100,	 /* maximum number of LSH iteration */
	BLK_SIZE = 2000, /* block size of DP-means  */
	TOT_SIZE = 0,	 /* number of sequence */
	SA_THRE = 3000,	 /* filter threshold of suffix array  */
	DP_ITER = 5,	 /* default iteration number of DP-means */
	THREAD_NUM = 2,  /* number of pthread */ 
	PARTITION = 1,   /* 1 for LSH partition, 2 for PCA partition */
	BIG_KMER = 32,   /* default value for BigKmer mapping */
	QSORT_BY = 1,
	KMER_ACC = 3,
	LOG_MODE = 0;

double	ID_THRE = 0.97,	   /*  clustering threshold  */
	LEVEL_THRE = 0.01; /*  convergence threshold for LSH iteration */


MPI_Status status;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER; 
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER; 

cls_log *logger = NULL;
vector<string> raw_seq;
vector<long> raw_hash;
long *seq = NULL;
int *seq_weight = NULL;
int seqlen;
double *pca_x = NULL, *pca_y = NULL;
double *pca_U = NULL;
int** pca_buf = NULL;

vector<list<int>* > blist;
vector<vector<vector<long> > > result;
vector<int> result_weight;
vector<boost::unordered_set<pair<long, long> > > klink;
boost::unordered_map<long, int> cl_id;
vector<vector<long> > cl;

int *displs = NULL, *rcounts = NULL;
long *seqbuf = NULL;

void print_usage(char* prog)
{
	if (MPI_RANK == 0) {
		printf("Usage: mpiexec -np [Host number] %s [Options]\n", prog);
		printf("\n");
		printf("  Note: Host number must be greater than 1: one for master, the others for slaves.\n");
		printf("\n");
		printf("Options\n");
		printf("   -i  input file in fasta format [required] \n");
		printf("   -o  output prefix [required]\n");
		printf("   -p  clustering threshold, default 0.97\n");
		printf("   -c  thread number, default 2\n");
		printf("   -b  block size of DP-means, default 2000\n");
		printf("   -l  maximum LSH iteration number, default 100\n");
		printf("   -x  convergence threshold for LSH iteration, default 0.01\n");
		printf("   -d  iteration number of DP-mean algorithm, default 5\n");
		printf("   -s  suffix array filter threshold, default 3000\n");
		printf("   -q  0/1, sort by length/weight in Big-Kmer Mapping, default 1\n");
		printf("   -h  print this help\n");
		printf("\n");
		printf("Example:\n");
		printf("   $ mpiexec -np 2 %s -i db.fa -o db_output -c 8\n", prog);
	}
	
	MPI::Finalize();
	exit(0);
}

long BKDRHash(const char *str)
{
	unsigned long seed = 1313; 
	unsigned long hash = 0;
 
	while (*str)
	{
		hash = hash * seed + (*str++);
	}
 
	return (long) hash;
}

void calc_kmer(int* buf, const string& str) 
{
	if (str.size() <= PCA_KMER) {
		logger->info("error, str.size() <= PCA_KMER");
		exit(-1);
	}
	memset(buf, 0, PCA_U_LEN * sizeof(int));

	int k = 0;
	for (int i = 0; i < PCA_KMER - 1; ++i) {
		switch(str[i]){
			case 'A': case 'a': k = k * 4 + 0; break;
			case 'T': case 't': k = k * 4 + 1; break;
			case 'C': case 'c': k = k * 4 + 2; break;
			case 'G': case 'g': k = k * 4 + 3; break;
			default:
				;//logger->info("[%c] warning not ATCG", str[i]); 
		}
	}
	
	for (int i = PCA_KMER - 1; i < str.size(); ++i) {
		switch(str[i]){
			case 'A': case 'a': k = k * 4 + 0; break;
			case 'T': case 't': k = k * 4 + 1; break;
			case 'C': case 'c': k = k * 4 + 2; break;
			case 'G': case 'g': k = k * 4 + 3; break;
			default:
				;//logger->info("[%c] warning not ATCG", str[i]); 
		}
		++buf[k %= PCA_U_LEN];
	}
	
	/*
	for (int i = 0; i < PCA_U_LEN; i++)
		printf("%d, ", buf[i]);
	printf(" size=[%d]\n", str.size());
	*/
}

double rnd01() {
	union {
		double d;
		unsigned long u;
	} x;
	x.u = ((unsigned long)rand() << 21) | (rand() & 0x00030000U) |  0x3FF0000000000000ULL;
	return x.d - 1.0;
}

void gen_projection() {

#ifdef USING_PCA

	if (PARTITION == 0) {
		logger->info("partitioning using PCA");

		/* 3x3 matrix A
		 * 76 25 11
		 * 27 89 51
		 * 18 60 32
		 */

		static doublereal X[PCA_U_LEN * PCA_SELECT]; // = {76, 27, 18, 25, 89, 60, 11, 51, 32};
		static doublereal S[PCA_U_LEN * PCA_U_LEN];  // = {0,0,0,0,0,0,0,0,0};

		static int kmer_buf[PCA_U_LEN];
		logger->info("X=%ld S=%ld kmer_buf=%ld, [%ld, %ld, %ld]", &X[0], &S[0], &kmer_buf[0], X, S, kmer_buf);
		double ratio = (double)PCA_SELECT / seqlen;
		int cnt = 0;
		for (int i = 0; i < seqlen && cnt < PCA_SELECT; ++i) {
			if (ratio > rnd01()) {
				calc_kmer(kmer_buf, raw_seq[seq[i]]);
				for (int j = 0; j < PCA_U_LEN; ++j) {
					X[cnt * PCA_U_LEN + j] = (double) kmer_buf[j] / (double) raw_seq[seq[i]].size();
					//printf("X = %lf\n", X[cnt * PCA_U_LEN + j]);
				}
				cnt++;
			}
		}

		integer dim1 = PCA_U_LEN, dim2 = cnt;
		doublereal alpha = 1.0 / cnt, zero = 0;

		char charN = 'N', charT = 'T', charV = 'V';
		dgemm_(&charN, &charT, &dim1, &dim1, &dim2, &alpha, X, &dim1, X, &dim1, &zero, S, &dim1);

		integer info ;

		static doublereal wr[PCA_U_LEN], wi[PCA_U_LEN];
		static doublereal vr[PCA_U_LEN * PCA_U_LEN], vl[PCA_U_LEN * PCA_U_LEN];

		integer lwork = PCA_U_LEN * PCA_U_LEN * 4 ;
		static doublereal work[PCA_U_LEN * PCA_U_LEN * 4 ];

		dgeev_(&charN, &charV, &dim1, S, &dim1, wr, wi, vl , &dim1 , vr, &dim1, work, &lwork, &info);

		for (int i = 0; i < PCA_U_LEN -1 ; i++) {
			if (wi[i] > 1e-5) {
				logger->info("error, wi>0");
			break;
			}
			if (wr[i] < wr[i + 1]) {
				logger->info("error wr asc: %lf < %lf", wr[i], wr[i + 1]);
			break;
			}
		}

		for (int i = 0; i < PCA_U_LEN * 2; i++) {
			pca_U[i] = vr[i];
		}
	}else
#endif

	{
		logger->info("partitioning using LSH");
		
		for (int i = 0; i < PCA_U_LEN * 2; i += 2) {
			double u1 = rnd01(), u2 = rnd01();
			while (u1 == 0.0 || u1 == 1.0) u1 = rnd01();
			while (u2 == 0.0 || u2 == 1.0) u2 = rnd01();
			pca_U[i] = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
			pca_U[i + 1] = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
		}
	}
}

void* do_hashing(void* parm) 
{
	int id = (long) parm;
	for (size_t i = id; i < raw_seq.size(); i += THREAD_NUM) {
		raw_hash[i] = BKDRHash(raw_seq[i].c_str());
	}
	return NULL;
}

void* do_pca(void* parm) 
{
	int id = (long) parm;
	int* buf = pca_buf[id];

	for (int p = id; p < seqlen; p += THREAD_NUM) {
		calc_kmer(buf, raw_seq[seq[p]]);
		pca_x[p] = pca_y[p] = 0;
		for (int i = 0; i < PCA_U_LEN; ++i) {
			pca_x[p] += buf[i] * pca_U[i];
			pca_y[p] += buf[i] * pca_U[i + PCA_U_LEN];
		}
		pca_x[p] /= (double) raw_seq[seq[p]].size();
		pca_y[p] /= (double) raw_seq[seq[p]].size();
	}
	return NULL;
}

int reqblockid() 
{
	pthread_mutex_lock(&mutex); 
	int ret = 1;
	MPI_Send(&ret, 1, MPI_INT, 0, MPI_TAG_REQ_BLOCKID, MPI_COMM_WORLD);
	MPI_Recv(&ret, 1, MPI_INT, 0, MPI_TAG_ANS_BLOCKID, MPI_COMM_WORLD, &status);  
	pthread_mutex_unlock(&mutex); 
	return ret;
}

void qsort_byid(long l , long r) {
	long m = seq[(l+r)/2], i = l, j = r;
	while (i <= j) 
	{
		while (seq[i] < m) i++;
		while (seq[j] > m) j--;
		if (i <= j) {
			swap(seq[i], seq[j]);
			swap(seq_weight[i], seq_weight[j]);
			i++, j--;
		}
	}
	if (l < j) qsort_byid(l, j);
	if (i < r) qsort_byid(i, r);
}

void qsort_bylen(long l , long r) {
	long m = raw_seq[seq[(l+r)/2]].size(), i = l, j = r;
	while (i <= j) 
	{
		while (raw_seq[seq[i]].size() < m) i++;
		while (raw_seq[seq[j]].size() > m) j--;
		if (i <= j) {
			swap(seq[i], seq[j]);
			swap(seq_weight[i], seq_weight[j]);
			i++, j--;
		}
	}
	if (l < j) qsort_bylen(l, j);
	if (i < r) qsort_bylen(i, r);
}

void qsort_byweight(long l , long r) {
	int m = seq_weight[(l+r)/2];
	long i = l, j = r;
	while (i <= j) 
	{
		while (seq_weight[i] < m) i++;
		while (seq_weight[j] > m) j--;
		if (i <= j) {
			swap(seq[i], seq[j]);
			swap(seq_weight[i], seq_weight[j]);
			i++, j--;
		}
	}
	if (l < j) qsort_byweight(l, j);
	if (i < r) qsort_byweight(i, r);
}

void ansblockid() {
	int alive = (MPI_SIZE - 1), ret, blockid = 0;
	while (alive > 0) {
		MPI_Recv(&ret, 1, MPI_INT, MPI_ANY_SOURCE, MPI_TAG_REQ_BLOCKID, MPI_COMM_WORLD, &status);
		int rank = status.MPI_SOURCE;
		if (ret == 1) {
			MPI_Send(&blockid, 1, MPI_INT, rank, MPI_TAG_ANS_BLOCKID, MPI_COMM_WORLD);
			++blockid;
		}else{
			if (ret != 0) {
				logger->info("error: ansblockid got ret != 0");
			}
			--alive;
		}
	}
	logger->info("ansblockid (master process) finished");
}

void* do_binning(void* parm)
{
	int tid = (long) parm,
		hid = (MPI_RANK - 1) * THREAD_NUM + tid,
		tot = (MPI_SIZE - 1) * THREAD_NUM;
	
	//logger->info("binning: pthread id=[%ld], handler id= [%d]", (long)parm, hid);

	boost::unordered_map<cls_kmer, vector<pair<long, int> >*, cls_kmer::hash_value > bin;

	for (long i = 0, g = 1; i < seqlen; i++) 
	{
		if (i == int(0.1 * seqlen * g)) {
			logger->info("binning progress: handler_id= [%d], %d\%", hid, g * 10);
			g++;
		}
		string& s = raw_seq[seq[i]];
		cls_kmer kmer(BIG_KMER, s.c_str());
		for (int j = BIG_KMER; j <= s.size(); j++) 
		{
			if (kmer % tot == hid) {
				if (bin.find(kmer) == bin.end()) {
					vector<pair<long, int> >* new_item = new vector<pair<long, int> >();
					new_item->push_back(make_pair<long, int> (i, j - BIG_KMER));
					bin[kmer] = new_item;
				}else{
					if (bin[kmer]->back().first != i)
						bin[kmer]->push_back(make_pair<long, int> (i, j - BIG_KMER));
				}
			}
			kmer.push_back(s[j]);
		}
	}

	cls_alignment *align = new cls_alignment(raw_seq);

	/*
	KMER_ACC = 0, 1, 2, 3
	  0: check all pair
	  1: check pair until find conn
	  2: check all center
	  3: check center until find conn
	*/
	for (boost::unordered_map<cls_kmer, vector<pair<long, int> >* >::iterator it = bin.begin(); it != bin.end(); ++it) 
	{
		if (it->second->size() > 1 && it->second->size() < 500) {
			if (KMER_ACC <= 1) {
			//logger->info("handler = %d, bin size = %d", hid, it->second->size()); 
				for (vector<pair<long, int> >::iterator jt = it->second->begin(); jt != it->second->end(); ++jt) {
					for (vector<pair<long, int> >::iterator kt = jt;  kt != it->second->end(); ++kt) 
						if (kt->first != jt->first)
						//seq may have repeat region
						{
							double d = align->get_similarity(seq[jt->first], seq[kt->first], ID_THRE);
							if (d >= ID_THRE) {
								klink[tid].insert(make_pair(jt->first, kt->first));
								//logger->info("jt= (%ld, %d) kt= (%ld, %d), d= %lf", seq[jt->first], jt->second, seq[kt->first], kt->second, d);
								if (KMER_ACC == 1) break;
							}
						}
				}
			}else {
			
				vector<long> ctr_list;
				vector<pair<long, int> >::iterator jt = it->second->begin();
				ctr_list.push_back(jt->first);
				int pre_id = jt->first;
				for (++jt; jt != it->second->end(); ++jt) {
					if (jt->first == pre_id)
						continue;
					for (int k = 0; k < ctr_list.size(); k++) {
						double d = align->get_similarity(seq[ctr_list[k]], seq[jt->first], ID_THRE);
						if (d >= ID_THRE) {
							klink[tid].insert(make_pair(ctr_list[k], jt->first));
							//logger->info("jt= (%ld, %d) kt= (%ld, %d), d= %lf", seq[jt->first], jt->second, seq[kt->first], kt->second, d);
							if (KMER_ACC == 3) break;
						}
					}
					pre_id = jt->first;
				}
			}
		}
	}

	logger->info("link size= %d", klink[tid].size());
	//sort(klink[tid].begin(), klink[tid].end());

	delete align;
}


void* do_clustering(void* parm) 
{
	long pid = (long) parm;
	logger->info("pid= [%ld] start for clustering", pid);

	int blockid;
	while ((blockid = reqblockid()) < blist.size()) {
		vector<string> raw;
		vector<long> raw_weight;
		vector<long> id;
		for (list<int>::iterator it = blist[blockid]->begin(); it != blist[blockid]->end(); ++it) {
			raw.push_back(raw_seq[seq[*it]]);
			raw_weight.push_back(seq_weight[*it]);
			//logger->info("raw[%d], weight=%d", id.size(), raw_weight.back());
			id.push_back(seq[*it]);
		}

		logger->info("pid= [%ld] finsh one block, block size= [%d], cluster number= [%d]", pid, blist[blockid]->size(), raw.size());

		cls_cluster *handler = new cls_cluster(raw, raw_weight, DP_ITER, SA_THRE, ID_THRE);

		if (SA_THRE == 0) {
			SA_THRE = handler->estimate_sa_threshold();
		}
		handler->process();

		pthread_mutex_lock(&mutex2); 
		
		int p = result.back().size();
		vector<int> subresult;
		handler->export_result(result.back(), result_weight);
		for (int i = p; i < result.back().size(); i++) {
			for (int j = 0; j < result.back()[i].size(); j++) {
				result.back()[i][j] = id[result.back()[i][j]];
			}
		}
		pthread_mutex_unlock(&mutex2); 

		delete handler;

	}

	return NULL;
}

void read_fasta() 
{
	logger->info("reading input [%s]", INPUT_FILE);

	FILE *fp = fopen(INPUT_FILE, "r");

	if (fp == NULL) {
		logger->fatal("can not open file [%s]", INPUT_FILE);
		MPI::Finalize();
		exit(-1);
	}

	if (MPI_RANK == 0) {
		// Master node won't read input file
		fclose(fp);
		return;
	}
	
	if (TOT_SIZE != 0)
		raw_seq.reserve(TOT_SIZE);
	else
		raw_seq.reserve(1 << 20);

	char buf[MAX_SEQ_LEN];
	char *fg = fgets(buf, MAX_SEQ_LEN, fp);
	while (!feof(fp) && (raw_seq.size() < TOT_SIZE || TOT_SIZE == 0)) 
	{
		int len = 0;
		while (!feof(fp)) {
			fg = fgets(buf + len, MAX_SEQ_LEN, fp);
			if (buf[len] == '>')
				break;
			len += strlen(buf + len) - 1;
		}
		buf[len] = '\0';
		raw_seq.push_back(buf);
	}

	if (TOT_SIZE == 0) 
		TOT_SIZE = raw_seq.size();
	
	fclose(fp);

	//unique
	logger->info("original sequences number= %d", raw_seq.size());
	
	logger->info("removing duplicated sequences");

	pthread_t tid[THREAD_NUM];

	raw_hash.resize(raw_seq.size());

	for (int i = 0; i < THREAD_NUM; i++)
		pthread_create(&tid[i], NULL, do_hashing, (void*) (long) i);
	
	for (int i = 0; i < THREAD_NUM; i++)
		pthread_join(tid[i], NULL);

	multimap<long, size_t> unique_map;
	seqlen = 0;
	
	for (size_t i = 0; i < raw_seq.size(); i++) {
		bool find = false;
		if (unique_map.find(raw_hash[i]) != unique_map.end()) {
			pair<multimap<long, size_t>::iterator, multimap<long, size_t>::iterator> range = unique_map.equal_range(raw_hash[i]);
			for (multimap<long, size_t>::iterator iter = range.first; iter != range.second; ++iter) {
				if (raw_seq[iter->second] == raw_seq[i]) {
					find = true;
					raw_hash[i] = -((long)iter->second);
					++raw_hash[iter->second];
					break;
				}
			}
		}
		if (!find) {
			unique_map.insert(pair<long, size_t>(raw_hash[i], i));
			raw_hash[i] = 1;
			++seqlen;
		}
	}

	logger->info("unique sequences number= %d", seqlen);

	seq = new long[seqlen];
	seq_weight = new int [seqlen];

	for (size_t i = 0, j = 0; i < raw_seq.size(); i++)
		if (raw_hash[i] > 0) {
			seq[j] = i;
			seq_weight[j] = raw_hash[i];
			++j;
		}
}

int process() 
{
	pthread_t tid[THREAD_NUM];

	logger->info("start clustering");
	
	pca_U = new double[PCA_U_LEN * 2];

	int loop = 0;
	int pre_otus_num = INT_MAX;

	if (MPI_RANK == 0) {
		displs = new int[MPI_SIZE], rcounts = new int[MPI_SIZE];
		seq = seqbuf = NULL;

	}else{

		pca_x = new double[seqlen];
		pca_y = new double[seqlen];

		pca_buf = new int* [THREAD_NUM];
		for (int i = 0; i < THREAD_NUM; ++i)
			pca_buf[i] = new int[PCA_U_LEN];
	}

	while (loop++ < LEVEL) 
	{
		logger->info("LSH iteration, loop id= %d", loop);

		//generate random projection vector
		if (MPI_RANK == 1) {
			logger->info("generating projection vector");
			gen_projection();
		}

		MPI_Bcast(pca_U, PCA_U_LEN * 2, MPI_DOUBLE, 1, MPI_COMM_WORLD);

		//perform data partition
		if (MPI_RANK != 0) {
				
				logger->info("calculating projection of K-mer vector");
				for (int i = 0; i < THREAD_NUM; i++)
					pthread_create(&tid[i], NULL, do_pca, (void*) (long) i);
				
				for (int i = 0; i < THREAD_NUM; i++)
					pthread_join(tid[i], NULL);
				
				logger->info("performing data partition");

				double max_x=pca_x[0], 
					   min_x=pca_x[0], 
					   max_y=pca_y[0], 
					   min_y=pca_y[0];

				for (int i = 0; i < seqlen; i++)
				{
					max_x = max(max_x, pca_x[i]);
					min_x = min(min_x, pca_x[i]);
					max_y = max(max_y, pca_y[i]);
					min_y = min(min_y, pca_y[i]);
				}

				int bedge = max(3, (int)sqrt(seqlen / BLK_SIZE * 2)),
					bcount = bedge * bedge; //bnum : edge 

				double delta_x = (max_x - min_x) / bedge;
				double delta_y = (max_y - min_y) / bedge;

				blist.clear();

				for (int i = 0; i < bcount; i++)
					blist.push_back(new list<int>());

				for (int i = 0; i < seqlen; i++) {
					
					int dx = (pca_x[i] - min_x) / delta_x;
					int dy = (pca_y[i] - min_y) / delta_y;

					if (dx == bedge) 
						--dx;
					if (dy == bedge) 
						--dy;

					int id = dx * bedge + dy;
					if (blist[id]->size() >= BLK_SIZE) {
						blist.push_back(blist[id]);
						blist[id] = new list<int>();
					}
					blist[id]->push_back(i);
				}

				int blast = 0;
				while (blast < blist.size()) {
					while (blist[blast]->size() >= BLK_SIZE && blast < blist.size()) 
						blast++;
					int i = blast + 1;
					while (blist[blast]->size() < BLK_SIZE && i < blist.size())
						blist[blast]->splice(blist[blast]->end(), *blist[i++]);
					if (i == blist.size()) 
						break;
					++blast;
				}

				while (blist.back()->size() == 0) {
					delete blist.back();
					blist.pop_back();
				}

		}


		//perform DP-means for clustering
		result.push_back(vector<vector<long> >());

		if (MPI_RANK == 0) {
			ansblockid();
			seqlen = 0;

		}else{

			logger->info("performing DP-means for clustering ...");
			for (int i = 0; i < THREAD_NUM; i++)
				pthread_create(&tid[i], NULL, do_clustering, (void*) (long) i);

			for (int i = 0; i < THREAD_NUM; i++)
				pthread_join(tid[i], NULL);
			
			int ret = 0;
			MPI_Send(&ret, 1, MPI_INT, 0, MPI_TAG_REQ_BLOCKID, MPI_COMM_WORLD);

			for (int i = 0; i < result.back().size(); i++) {
				//centroid were used in next LSH iteration
				seq[i] = (result.back())[i].front();
				seq_weight[i] = result_weight[i];
			}
			seqlen = result.back().size();
		}


		//sync and merge all clusters
		logger->info("sync clustering result of each block");

		MPI_Gather(&seqlen, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//logger->info("MPI_Gather seqlen finished");

		int buflen = 0;

		if (MPI_RANK == 0) {
			displs[0] = 0;
			for (int i = 1; i < MPI_SIZE; ++i) {
				displs[i] = displs[i-1] + rcounts[i-1];
			}
			
			buflen = displs[MPI_SIZE - 1] + rcounts[MPI_SIZE - 1];
			
			if (seq == NULL) {
				seq = new long [buflen];
				seq_weight = new int [buflen];
				seqbuf = new long[buflen];
			}

		}

		MPI_Gatherv(seq, seqlen, MPI_LONG, seqbuf, rcounts, displs, MPI_LONG, 0, MPI_COMM_WORLD);
		//logger->info("MPI_Gatherv seq finished");

		if (MPI_RANK == 0) {
			memcpy(seq, seqbuf, sizeof(long) * buflen);
		}
		MPI_Gatherv(seq_weight, seqlen, MPI_INT, seqbuf, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
		//logger->info("MPI_Gather seq_weight finished");

		if (MPI_RANK == 0) {
			memcpy(seq_weight, seqbuf, sizeof(int) * buflen);

			seqlen = buflen;
			//for FIX ANS, sort seq[] and seq_weight[]
			qsort_byid(0, seqlen - 1);
		}

		MPI_Bcast(&seqlen, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//logger->info("MPI_Bcast seqlen finished");

		MPI_Bcast(seq, seqlen, MPI_LONG, 0, MPI_COMM_WORLD);
		//logger->info("MPI_Bcast seq finished");

		MPI_Bcast(seq_weight, seqlen, MPI_INT, 0, MPI_COMM_WORLD);
		//logger->info("MPI_Bcast seq_weight finished");

		int singleton = 0;
		for (int i = 0; i < seqlen; i++)
			singleton += (seq_weight[i] == 1);

		logger->info("LSH iteration finish, cluster num= [%d], singleton= [%d],  nonsingleton= [%d]", seqlen, singleton, seqlen - singleton);

		//clearn up 
		//logger->info("start clean up");
		for (int i = 0; i < blist.size(); i++)
			delete blist[i];

		blist.clear();
		result_weight.clear();
		
		double percent = (double) (pre_otus_num - seqlen) / (double) seqlen;
		if (percent < LEVEL_THRE) {
			logger->info("previous cluster num= %d, merge percentage= %.4f, LSH approach convergence", pre_otus_num, percent);
			break;

		}else{
			pre_otus_num = seqlen;
		}
	}

	//free memory 
	delete [] pca_U;

	if (MPI_RANK == 0) {
		delete [] displs;
		delete [] rcounts;
		displs = rcounts = NULL;

		if (seq != NULL) {
			delete [] seq;
			delete [] seq_weight;
			delete [] seqbuf;
			seq = NULL;
			seq_weight = NULL;
			seqbuf = NULL;
		}
	}else{
		delete [] pca_x;
		delete [] pca_y;
		pca_x = pca_y = NULL;

		for (int i = 0; i < THREAD_NUM; ++i)
			delete [] pca_buf[i];
		delete [] pca_buf;
		pca_buf = NULL;
	}

	//Big-Kmer Mapping
	if (KMER_ACC < 4) {

		logger->info("start Big-Kmer Mapping clustering, seqlen= %d", seqlen);

		int send_len = 0, recv_len = 0;
		long *send_buf = NULL, *recv_buf = NULL;

		if (MPI_RANK != 0) {

			if (QSORT_BY == 0)
				qsort_bylen(0, seqlen - 1);
			else
				qsort_byweight(0, seqlen - 1);

			klink.resize(THREAD_NUM);

			for (int i = 0; i < THREAD_NUM; i++)
				pthread_create(&tid[i], NULL, do_binning, (void*) (long) i);
			
			for (int i = 0; i < THREAD_NUM; i++)
				pthread_join(tid[i], NULL);
			
			//sync results

			logger->info("start sync Big-Kmer Mapping result");

			for (int i = 0; i < THREAD_NUM; ++i)
				send_len += klink[i].size();

			send_buf = new long [send_len * 2];
			
			send_len = 0;
			for (int i = 0; i < THREAD_NUM; ++i)
				for (boost::unordered_set<pair<long, long> >::iterator jt = klink[i].begin(); jt != klink[i].end(); ++jt) 
				{
					send_buf[send_len++] = jt->first;
					send_buf[send_len++] = jt->second;
				}

		}

		MPI_Barrier(MPI_COMM_WORLD);

		if (MPI_RANK == 1) {
			displs = new int[MPI_SIZE], rcounts = new int[MPI_SIZE];
		}

		//logger->info("send_len= %d", send_len);
		MPI_Gather(&send_len, 1, MPI_INT, rcounts, 1, MPI_INT, 1, MPI_COMM_WORLD);

		if (MPI_RANK == 1) {

			displs[0] = 0;
			for (int i = 1; i < MPI_SIZE; ++i) {
				displs[i] = displs[i-1] + rcounts[i-1];
			}
			
			recv_len = displs[MPI_SIZE - 1] + rcounts[MPI_SIZE - 1];
			//logger->info("recv_len= %d", recv_len);
			recv_buf = new long[recv_len];
		}

		logger->info("Gathering Big-Kmer Mapping result");
		MPI_Gatherv(send_buf, send_len, MPI_LONG, recv_buf, rcounts, displs, MPI_LONG, 1, MPI_COMM_WORLD);

		if (MPI_RANK == 1) 
		{
			vector<pair<long, long> > tlink;
			for (int i = 0; i < recv_len; i += 2) {
				tlink.push_back(make_pair(recv_buf[i], recv_buf[i + 1]));
			}
			sort(tlink.begin(), tlink.end());
			
			for (size_t i = 0, j = 0; i < seqlen; i++) {
				if (cl_id.find(seq[i]) == cl_id.end()) {
					cl_id[seq[i]] = cl.size();
					cl.push_back(vector<long>());
					cl.back().push_back(seq[i]);
				}
				long now_id = cl_id[seq[i]];

				while (j < tlink.size() && tlink[j].first == i) {
					long t = tlink[j].second;
					if (cl_id.find(seq[t]) == cl_id.end()) {
						cl_id[seq[t]] = now_id;
						cl[now_id].push_back(seq[t]);
					}
					++j;
				}
			}

			/*
			for (size_t i = 0; i < raw_seq.size(); i++)
				if (raw_hash[i] <= 0) {
					cl[cl_id[-raw_hash[i]]].push_back(i);
				}

			fp = fopen(OUTPUT_FILE, "w");
			for (size_t i = 0; i < cl.size(); i++) {
				sort(cl[i].begin(), cl[i].end());
				for (size_t j = 0; j < cl[i].size(); j++)
					fprintf(fp, "%ld ", cl[i][j]);
				fprintf(fp, "\n");
			}
			fclose(fp);
			*/
		}

		logger->info("Big-Kmer Mapping finished");

		//clean up
		if (send_buf != NULL) delete [] send_buf;
		if (recv_buf != NULL) delete [] recv_buf;
	}

	//construct final answer
	logger->info("start constructing final answer");

	if (result.size() > 0) {

		if (MPI_RANK != 0) {
			if (MPI_RANK == 1) {
				rcounts = new int[MPI_SIZE], displs = new int [MPI_SIZE];
				seqbuf = new long [(long)(TOT_SIZE * 1.5)];
			}	
			delete [] seq;
			delete [] seq_weight;

			int l = result[0].size();
			for (int k = 0; k < result[0].size(); k++) 
				l += result[0][k].size();
			seq = new long [l + 1];
		}

		for (int lv = result.size() - 1; lv >= 0; --lv) {
			seqlen = 0;
			for (int k = 0; k < result[lv].size(); k++) {
				for (int i = 0; i < result[lv][k].size(); i++) {
					seq[seqlen++] = result[lv][k][i];
				}
				seq[seqlen++] = -1;
			}

			MPI_Gather(&seqlen, 1, MPI_INT, rcounts, 1, MPI_INT, 1, MPI_COMM_WORLD);
			
			if (MPI_RANK == 1) {
				displs[0] = 0;
				for (int i = 1; i < MPI_SIZE; ++i) {
					displs[i] = displs[i-1] + rcounts[i-1];
				}
			}

			MPI_Gatherv(seq, seqlen, MPI_LONG, seqbuf, rcounts, displs, MPI_LONG, 1, MPI_COMM_WORLD);

			if (MPI_RANK == 1) {
				int l = displs[MPI_SIZE - 1] + rcounts[MPI_SIZE - 1];
				for (int i = 0; i < l; i++) {
					
					if (lv == result.size() - 1 && KMER_ACC > 3) {
						cl_id[seqbuf[i]] = cl.size();
						cl.push_back(vector<long>());
						cl.back().push_back(seqbuf[i]);
					}
					int cid = cl_id[seqbuf[i]];
					while (seqbuf[++i] != -1) {
						cl[cid].push_back(seqbuf[i]);
						cl_id[seqbuf[i]] = cid;
					}
				}
			}
			/*
			if (MPI_RANK == 1) {
				for (int i = 0; i < cl.size(); i++) {
					for (int j = 0; j < cl[i].size(); j++)
						printf("%d ", cl[i][j]);
					printf("\n");
				}
				printf("=============================\n");	
			}
			*/
		}

		if (MPI_RANK != 0) {
			if (MPI_RANK == 1) {
				delete [] rcounts;
				delete [] displs;
				delete [] seqbuf;
			}	
			delete [] seq;
		}

	}

	long count[2];

	if (MPI_RANK == 1) {
		for (size_t i = 0; i < raw_seq.size(); i++) {
			if (raw_hash[i] <= 0) {
				cl[cl_id[-raw_hash[i]]].push_back(i);
			}
		}

		string clstr = string(OUTPUT_FILE) + ".clstr";
		fp = fopen(clstr.c_str(), "w");
		for (size_t i = 0; i < cl.size(); i++) {
			fprintf(fp, ">%ld\n%s\n", cl[i][0], raw_seq[cl[i][0]].c_str());
		}
		fclose(fp);	

		count[0] = cl.size();
		count[1] = 0;

		string fpout = string(OUTPUT_FILE) + ".ans";
		fp = fopen(fpout.c_str(), "w");
		for (size_t i = 0; i < cl.size(); i++) {
			sort(cl[i].begin(), cl[i].end());
			for (size_t j = 0; j < cl[i].size(); j++)
				fprintf(fp, "%ld ", cl[i][j]);
			count[1] += (cl[i].size() == 1);
			fprintf(fp, "\n");
		}
		fclose(fp);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(count, 2, MPI_LONG, 1, MPI_COMM_WORLD);
	if (MPI_RANK == 0) {
		logger->info("final cluster num= [%ld], singleton= [%ld], nonsingleton= [%ld]", count[0], count[1], count[0] - count[1]);
	}

}

int main(int argc, char **argv)
{

	clock_t start = clock();

	MPI_Init(&argc, &argv);

	logger = cls_log::get_instance();
	
	//MPI_SIZE = MPI::COMM_WORLD.Get_size();
	MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
	MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
	//MPI_RANK = MPI::COMM_WORLD.Get_rank();

#ifdef FIX_RAND_SEED
	srand(1031);
#else
	srand(time(NULL));
#endif
		
	int ch;
	while ((ch = getopt(argc, argv, "i:o:b:t:l:x:d:p:s:c:r:k:q:a:u:h")) != -1) {  
		switch (ch) { 
			case 'i': INPUT_FILE = optarg; break; 
			case 'o': OUTPUT_FILE = optarg; break; 
			case 'b': BLK_SIZE = atoi(optarg); break; 
			case 't': TOT_SIZE = atoi(optarg); break; 
			case 'l': LEVEL = atoi(optarg); break; 
			case 'x': LEVEL_THRE = atof(optarg);break;
			case 'd': DP_ITER = atoi(optarg); break;
			case 'p': ID_THRE = atof(optarg); break;
			case 's': SA_THRE = atoi(optarg); break;
			case 'c': THREAD_NUM = atoi(optarg); break;
			case 'r': PARTITION = atoi(optarg); break;
			case 'k': BIG_KMER = atoi(optarg); break;
			case 'q': QSORT_BY = atoi(optarg); break;
			case 'a': KMER_ACC = atoi(optarg); break;
			case 'u': LOG_MODE = atoi(optarg);  break;
			case 'h': print_usage(argv[0]);  break;
			default:
				break;
		}  
	}

	if (MPI_SIZE == 1) {
		logger->fatal("There should be at least 2 MPI processes. Use '%s -h' for more details", argv[0]);
		exit(-1);
	}

	if (!INPUT_FILE || !OUTPUT_FILE) {
		if (MPI_RANK == 0) {
			print_usage(argv[0]);
		}
	}

	logger->set_mode(LOG_MODE);

	string running_parameter = "";
	for (int i = 0; i < argc; i++)
		running_parameter += string(argv[i]) + " ";
	logger->info("Running parameter: %s", running_parameter.c_str());

	//logger->info("Rank[%d] startup", rank);

	read_fasta();
	process();

	logger->info("MPI_Rank[%d] return", MPI_RANK);
	
	if (MPI_RANK == 0) {
		logger->info("Timecost= [%.2lf]s", COUNT_SECONDS(start));
	}

	MPI::Finalize();
	return 0;
}
