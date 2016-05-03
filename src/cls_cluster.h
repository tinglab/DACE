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

#ifndef __cls__cluster__H__
#define __cls__cluster__H__

#include <string>
#include <vector>
#include <list>

#include "cls_alignment.h"
#include "cls_suffix_array.h"

//#define NO_SUFFIX_ARRAY
//#define SINGLE_ITER
#define LONGEST_CENTER

#define ERASE_CENTER
//#define RAND_SAMPLE_CENTER

using std::vector;

class cls_cluster
{
public:

	class st_center;

	class st_elem {
	public:
		st_elem(int id, int weight);
		~st_elem();
		void union_elem(st_elem* child);		

	public:
		enum elem_status {
			ACTIVE = 0,
			FIXED,
			UNIONED
		} status;

		int n, sa_begin, sa_end, weight;
		std::list<int> union_list;
		st_center* ctr;
		double similar;
	};

	class st_center {
	public:
		st_center(st_elem* center);
		~st_center();
		bool recenter(cls_alignment& align, double SPECIES_CALC);

	public:
		int hash, size;
		st_elem* one;
		std::vector<st_elem*> grp;

	};

	cls_cluster(std::vector<std::string>& raw, std::vector<long>& raw_weight, int iter, int sa_threshold, double id_threshold);

	~cls_cluster();

	void scan_center(st_center* c, std::list<st_elem*>::iterator scan_begin, std::list<st_elem*>::iterator scan_end);

	void process();

	void print_result();

	void export_result(vector<vector<long> >& result, vector<int>& weight) ;

	int estimate_sa_threshold();

private:
	cls_suffix_array sa;
	cls_alignment align;
	std::list<st_center*> center_list, fix_list;
	std::list<st_elem*> elem_list;
	int *hx, *hy;
	
	std::vector<std::string> raw;

	int ITER, SA_THRE;

	cls_log* logger;

public:
	double UNION_THRE, SPECIES_SOFT, SPECIES_HARD, SPECIES_CALC;
};

#endif
