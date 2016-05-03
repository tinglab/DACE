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
 
#include "cls_suffix_array.h"
#include "cls_common.h"

#include <iostream>
#include <cstdio>

using namespace std;

cls_suffix_array::cls_suffix_array()
{
	is_build = false;
	sa = height = wx = wy = bar = nrank = NULL;
}

cls_suffix_array::~cls_suffix_array()
{
	if (bar != NULL) {
		delete [] bar;
		bar = NULL;
	}
	if (sa != NULL) {
		delete [] sa;
		sa = NULL;
	}
	if (wx != NULL) {
		delete [] wx;
		wx = NULL;
	}
	if (wy != NULL) {
		delete [] wy;
		wy = NULL;
	}
	if (nrank != NULL) {
		delete [] nrank;
		nrank = NULL;
	}
	if (height != NULL) {
		delete [] height;
		height = NULL;
	}
}


int cls_suffix_array::size()
{
	return data.size();
}


int* cls_suffix_array::get_height()
{
	return is_build ? height : NULL;
}


int* cls_suffix_array::get_rank()
{
	return is_build ? wx : NULL;
}


int* cls_suffix_array::get_sa()
{
	return is_build ? sa : NULL;
}

int cls_suffix_array::cmp(int in[], int a, int b, int l)
{
	return in[a] == in[b] && in[a + l] == in[b + l];
}

int cls_suffix_array::add(const string& str, long ptr)
{
	for (int i = 0; i < str.size(); i++) {
		switch(str[i]) {
			case 'A': case 'a': data.push_back(1); break;
			case 'T': case 't': data.push_back(2); break;
			case 'C': case 'c': data.push_back(3); break;
			case 'G': case 'g': data.push_back(4); break;
			default:
				data.push_back(1);
				//data.push_back(rand() % 4);
				//log("ERROR: rna=[%s] rna[i]=[%c]", str.c_str(), str[i]);
				//return -1;
		}
	}
	data.push_back(0);

	//binsearch can be used to save memory
	data_ptr.insert(data_ptr.end(), data.size() - data_ptr.size(), ptr);
	return 0;
}

long cls_suffix_array::find_item(int index)
{
	return index >= 0 && index < data_ptr.size() ? data_ptr[index] : -1;
}

int cls_suffix_array::build()
{
	int len = data.size();
	int m = 5;

	bar = new int [len];
	sa = new int [len];
	wx = new int [len];
	wy = new int [len];
	nrank = new int [len];
	height = new int [len];

	int i, *rank = wx, *sa_y = wy, *t;

	for (i = 0; i <= m; i++)
		bar[i] = 0;

	for (i = 0; i < len; i++)
		bar[rank[i] = data[i]]++;
	
	for (i = 0; i < m; i++)
		bar[i + 1] += bar[i];
	
	for (i = len-1; i >= 0; i--)
		sa[--bar[rank[i]]] = i;
	
	for (int p = 1, k = 1; p < len; k *= 2, m = p)
	{
		for (p = 0, i = len - k; i < len; i++)
			sa_y[p++] = i;

		for (i = 0; i < len; i++)
			if (sa[i] >= k)
				sa_y[p++] = sa[i] - k;

		for (i = 0; i < len; i++)
			nrank[i] = rank[sa_y[i]];

		for (i = 0; i <= m; i++)
			bar[i] = 0;

		for (i = 0; i < len; i++)
			bar[nrank[i]] ++;

		for (i = 0; i < m; i++)
			bar[i+1] += bar[i];

		for (i = len-1; i >= 0; i--)
			sa[--bar[nrank[i]]] = sa_y[i];

		for (t = sa_y, sa_y = rank, rank = t, p = 1, rank[sa[0]] = 0, i = 1; i < len; i++)
			rank[sa[i]] = cmp(sa_y, sa[i-1], sa[i], k) ? p-1 : p++;
	}

	

	for (i = 1; i < len; i++)
		rank[sa[i]] = i;

	int j = 0, k = 0;
	for (i = 0; i < len-1; i++) {
		for (k ? k-- : 0, j = sa[rank[i] - 1]; data[i + k] == data[j + k]; k++) ;
		height[rank[i]] = k;
	}

	is_build = true;
	return 0;
}
