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

#include "cls_kmer.h"
#include <cstring>
#include <cstdio>

cls_kmer::cls_kmer(int length, const char* s) 
{
	//assert(length > 0);
	kmerlen = length;
	size = (length - 1)/ 32 + 1;
	array = new unsigned long [size];
	memset(array, 0, sizeof(unsigned long) * size);

	for (int i = 0, p, q; i < length && *s; i++) {
		p = i / 32;
		q = (i % 32) * 2;
		switch(*s++){
			case 'A': case 'a': array[p] |= (0L << q); break;
			case 'T': case 't': array[p] |= (1L << q); break;
			case 'C': case 'c': array[p] |= (2L << q); break;
			case 'G': case 'g': array[p] |= (3L << q); break;
			default:
				;//log("[%c] warning not ATCG", str[i]); 
		}
	}
}

cls_kmer::~cls_kmer()
{
	delete [] array;
	array = NULL;
}

void cls_kmer::push_back(char c)
{
	unsigned long m = 0;
	switch(c){
		case 'A': case 'a': m = 0; break;
		case 'T': case 't': m = 1; break;
		case 'C': case 'c': m = 2; break;
		case 'G': case 'g': m = 3; break;
		default:
			;//log("[%c] warning not ATCG", str[i]); 
	}
	for (int i = size - 1, m_next, q = ((kmerlen - 1) % 32) * 2; i >= 0; --i, q = 62) {
		m_next = array[i] & 3;
		array[i] >>= 2;
		array[i] |= (m << q);
		m = m_next;
	}
}

void cls_kmer::print() const
{
	for (int i = 0, p, q, m; i < kmerlen; i++) {
		p = i / 32;
		q = (i % 32) * 2;
		m = (array[p] >> q) & 3;
		switch(m) {
			case 0: printf("A"); break;
			case 1: printf("T"); break;
			case 2: printf("C"); break;
			case 3: printf("G"); break;
			default: ;
		}
	}
	printf("\n");
}

bool cls_kmer::operator< (const cls_kmer &b) const
{
	//assert(size == b.size);
	for (int i = 0; i < size; i++)
		if (array[i] < b.array[i])
			return true;
		else if (array[i] > b.array[i])
			return false;
	return false;
}

long cls_kmer::operator% (const long m) const
{
	return array[0] % m;
}

bool cls_kmer::operator== (const cls_kmer &b) const
{
	//assert(size == b.size);
	for (int i = 0; i < size; i++)
		if (array[i] != b.array[i])
			return false;
	return true;
}


size_t cls_kmer::hash_value::operator() (const cls_kmer & a) const
{
	size_t h = 0;
	for (int i = 0; i < a.size; i++)
		h ^= a.array[i];
	return h;
}

cls_kmer::cls_kmer(const cls_kmer &b)
{
	kmerlen = b.kmerlen;
	size = b.size;
	array = new unsigned long [size];
	memcpy(array, b.array, sizeof(unsigned long) * size);
}


/*
int main()
{
	string a = "GTGCGCCAAACCGGAAAGACGGATGTATATGACCTGAGTATCTTAATGGGACATGAAAGCATCGAGACAACAAAAAAATATGAGCATTTCGCCCATGAAATCATTGCAGCAGAAAACAATATTTCACACTTGGATATGTGTCTAGGAAATGGACAAAAAACAGGGCGCGGGTAG";
	int k = 80;
	cls_kmer t(k, a.c_str());
	t.print();
	for (int i = k; i < a.size(); i++) {
		t.push_back(a[i]);
		t.print();
	}
}
*/