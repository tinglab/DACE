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

#ifndef __cls__suffix__array__H__
#define __cls__suffix__array__H__

#include <string>
#include <vector>

class cls_suffix_array
{
public:
	cls_suffix_array();
	~cls_suffix_array();

	int size();
	int add(const std::string& str, long ptr = 0);
	int build();

	int* get_height();
	int* get_rank();
	int* get_sa();

	long find_item(int index);

private:

	int cmp(int in[], int a, int b, int l);

	bool is_build;

	int *sa, *height, *wx, *wy, *bar, *nrank;

	std::vector<int> data;
	std::vector<long> data_ptr;
};

#endif
