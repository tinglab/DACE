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
#include <unistd.h>

#include "cls_common.h"
#include "cls_log.h"
#include "cls_alignment.h"
#include "cls_suffix_array.h"
#include "cls_cluster.h"

using namespace std;

cls_cluster::st_elem::st_elem(int id, int _weight)
	: n(id), weight(_weight), ctr(NULL), similar(0), status(ACTIVE)
{
}

cls_cluster::st_elem::~st_elem()
{
}
	
void cls_cluster::st_elem::union_elem(st_elem* child)
{
	child->status = UNIONED;
	this->weight += child->weight;
	this->union_list.push_back(child->n);
	this->union_list.splice(this->union_list.end(), child->union_list);
}


cls_cluster::st_center::st_center(st_elem* center)
	: one(center), hash(0xFFFF), size(center->weight)
{
}

cls_cluster::st_center::~st_center()
{
}

#ifndef RAND_SAMPLE_CENTER
bool grp_cmp(const cls_cluster::st_elem* a, const cls_cluster::st_elem* b) 
{
	return a->weight > b->weight;
}
#endif

bool cls_cluster::st_center::recenter(cls_alignment& align, double SPECIES_CALC)
{
	double max_sum = 0, min_sum = -1;
	double p = min(1.0, max(10.0 / (double)grp.size(), 0.1)); //select 10% (or at least 10) seqs
	int tmp_hash = 0;
	int calc_n=0, old_ctr=one->n;
	size = 0;
#ifndef RAND_SAMPLE_CENTER
	sort(grp.begin(), grp.end(), grp_cmp);
	int pnum = (int)(p * (double)grp.size());
	//printf("p=[%lf], pnum=[%d]", p, pnum);
	for (vector<st_elem*>::iterator it = grp.begin(); it != grp.end() && pnum; it++, --pnum) {
  	{
#else
	for (vector<st_elem*>::iterator it = grp.begin(); it != grp.end(); it++) {
		//calc recenter:
		if (rand01() < 1.0 - pow(1.0-p, (*it)->weight)) {
#endif
			calc_n++;
			double sum = 0;
			for (vector<st_elem*>::iterator jt = grp.begin(); jt != grp.end(); jt++) {
				if ((*jt)->n != (*it)->n){
					double dist = 1 - align.get_similarity((*it)->n, (*jt)->n, SPECIES_CALC);
					sum += dist * ((*jt)->weight);
				}
			}
			if (min_sum < 0 || sum < min_sum || sum == min_sum && one->weight < (*it)->weight) {
				one = (*it);
				min_sum = sum;
			}
			/*
			//logger->info("max_sum[%d]=%lf, new_sum[%d]=%lf", one->n, max_sum, (*it)->n, sum);			
			//if (one->n == (*it)->n) sum = 0;
			if (sum > max_sum || sum == max_sum && one->weight < (*it)->weight) {
				//printf("old_sum=%.5lf, new_sum=%.5lf, old->w=[%d], new->w=[%d], size=[%d]\n", max_sum, sum, one->weight, (*it)->weight, (int)grp.size());
				one = (*it);
				max_sum = sum;
			}
			*/
		}
		//calc hash:
		tmp_hash = ((tmp_hash + (*it)->n * 131171) ^ ((*it)->n * 79991));// + rand();

		//calc size
		size += (*it)->weight;
	}

	//logger->info("recenter, old_ctr=[%d] new_ctr=[%d] sample_n=[%d] grpsize=[%d] max_sum=[%.2lf]", old_ctr, one->n, calc_n, grp.size(), max_sum);

	tmp_hash ^= one->n;
	if (tmp_hash == hash) {
		for (vector<st_elem*>::iterator it = grp.begin(); it != grp.end(); it++) {
			(*it)->status = st_elem::FIXED;
		}
		return true;
	}else{
		one->ctr = this;
		hash = tmp_hash;
		return false;
	}
}


cls_cluster::cls_cluster(vector<string>& _raw, vector<long>& raw_weight, int iter, int sa_threshold, double id_threshold) 
	:raw(_raw), align(raw), ITER(iter), SA_THRE(sa_threshold)
{

	logger = cls_log::get_instance();

	SPECIES_SOFT = id_threshold;
	SPECIES_HARD = 1.0 - (1.0 - id_threshold) / 2.0;
	SPECIES_CALC = 1.0 - (1.0 - id_threshold) * 2.0;
	//UNION_THRE = max(0.985, SPECIES_HARD);
	UNION_THRE = 0.98;

	for (int i = 0; i < raw.size(); i++) {
		elem_list.push_back(new st_elem(i, raw_weight[i]));
		
#ifndef NO_SUFFIX_ARRAY
		elem_list.back()->sa_begin = sa.size();
		sa.add(raw[i]);
		elem_list.back()->sa_end = sa.size();
#endif
	}
#ifndef NO_SUFFIX_ARRAY
	sa.build();
#endif

	hx = hy = NULL;
}

cls_cluster::~cls_cluster()
{
	if (hx != NULL) {
		delete [] hx;
		hx = NULL;
	}
	if (hy != NULL) {
		delete [] hy;
		hy = NULL;
	}
}


bool estimate_sa_threshold_cmp(const pair<double, int>& a, const pair<double, int>& b)
{
	return a.first > b.first || a.first == b.first && a.second > b.second;
}

int cls_cluster::estimate_sa_threshold()
{
#ifdef NO_SUFFIX_ARRAY	
	return -1;
#endif

	logger->info("estimating suffix array threshold...");

	int* rank = sa.get_rank();
	int* height = sa.get_height();

	hx = new int [sa.size() + 10];
	hy = new int [sa.size() + 10];

	memset(hx, 0, sizeof(int) * sa.size());
	memset(hy, 0, sizeof(int) * sa.size());

	vector<int> val,val_all;

	for (int i = elem_list.front()->sa_begin; i < elem_list.front()->sa_end; i++) {
		if (rank[i] + 1 < sa.size())
			hx[rank[i] + 1] = height[rank[i] + 1];
		if (rank[i] - 1 >= 0) 
			hy[rank[i] - 1] = height[rank[i] - 1];
	}

	for (int i = 1, j = sa.size() - 2; i < sa.size(); i++, j--) {
		if (hx[i] == 0) hx[i] = min(height[i], hx[i - 1]);
		if (hy[j] == 0) hy[j] = min(height[j], hy[j + 1]);
	}

	for (list<st_elem*>::iterator it = elem_list.begin(); it != elem_list.end(); ++it) {
		st_elem* p = *it;

		if (p->n == elem_list.front()->n) continue;

		int s = 0;
		for (int i = p->sa_begin; i < p->sa_end; i++) 
			s += max(hy[rank[i]], hx[rank[i]]);

		double d = align.get_similarity(elem_list.front()->n, p->n, SPECIES_SOFT);
		
		if (d >= SPECIES_SOFT) {
			val.push_back(s);
		}

		if (d != 0) val_all.push_back(s);
	}

	if (val.size() ==0 && val_all.size() == 0) {
		logger->info("Warning: too small block, suffix array cannot determine threshold, running with out suffix array");
		return 0;
	}

	int ret;
	if (val.size() > 0) {
		sort(val.begin(), val.end());
		ret = val[min(10, (int)(val.size() * 0.05))];
	}else{
		sort(val_all.begin(), val_all.end());
		logger->info("Warning: too small block, suffix array may not accurate");
		ret = val_all[val_all.size() * 0.95];
	}
	
	logger->info("Estimate SA_THRE=[%d]", ret);
	
	delete [] hx; 
	delete [] hy;

	hx = hy = NULL;

	return SA_THRE = ret;

}

void cls_cluster::scan_center(st_center* c, list<st_elem*>::iterator scan_begin, list<st_elem*>::iterator scan_end)
{
	//logger->info("scan center in");
	
	clock_t start = clock();
	c->grp.clear();
	c->grp.push_back(c->one);
	c->one->ctr = c;
	c->one->similar = 1;

#ifndef NO_SUFFIX_ARRAY
	int* rank = sa.get_rank();
	int* height = sa.get_height();

	memset(hx, 0, sizeof(int) * sa.size());
	memset(hy, 0, sizeof(int) * sa.size());
	
	for (int i = c->one->sa_begin; i < c->one->sa_end; i++) {
		if (rank[i] + 1 < sa.size())
			hx[rank[i] + 1] = height[rank[i] + 1];
		if (rank[i] - 1 >= 0) 
			hy[rank[i] - 1] = height[rank[i] - 1];
	}

	for (int i = 1, j = sa.size() - 2; i < sa.size(); i++, j--) {
		if (hx[i] == 0) hx[i] = min(height[i], hx[i - 1]);
		if (hy[j] == 0) hy[j] = min(height[j], hy[j + 1]);
	}
#endif

	for (; scan_begin != scan_end; scan_begin++) {
		st_elem* p = *scan_begin;

		if (p->ctr != NULL && p->similar >= SPECIES_HARD) {
			//dist_ignore++;
			//logger->info("dist_ignore, p->n=%d", p->n);
			continue;
		}
		else if (p->ctr == NULL || p->ctr->one != p) {
#ifndef NO_SUFFIX_ARRAY
			int s = 0;
			if (SA_THRE > 0)
				for (int i = p->sa_begin; i < p->sa_end; i++) 
					s += max(hy[rank[i]], hx[rank[i]]);

			if (s >= SA_THRE) {
#endif
				double d = align.get_similarity(c->one->n, p->n, SPECIES_SOFT);
				double bias = 0;
				if (c->size != p->weight)
					bias = abs(log(log((double)max(c->size, p->weight)/(double)min(c->size,p->weight)))) / 100;
				//d += bias;
				//logger->info("size=%d, weight=%d, bias=%.4f, pd=%.4f, d=%.4f", c->size, p->weight, bias, d - bias, d);
				if (d > p->similar) {
					p->ctr = c;
					p->similar = d;
				}
#ifndef NO_SUFFIX_ARRAY	
			}else{
				//nfilter++;
			}
#endif
		}
	}
	//logger->info("scan center out, timecost=[%.2lf]s", COUNT_SECONDS(start));
}


void cls_cluster::process()
{
#ifdef LONGEST_CENTER

	for (list<st_elem*>::iterator it = elem_list.begin(); it != elem_list.end(); ++it) {
		for (list<st_elem*>::iterator jt = it; jt != elem_list.end(); ++jt) {
			if (raw[(*jt)->n].size() > raw[(*it)->n].size()) {
				swap(*jt, *it);
			}
		}
	}

#endif

	center_list.push_back(new st_center(elem_list.front()));

	hx = new int [sa.size() + 10];
	hy = new int [sa.size() + 10];

	int iter = ITER;

#ifdef SINGLE_ITER
	iter = 1;
#endif

	while (iter--) {
		clock_t start = clock();
		logger->info("-------------------ITER=%d Start!-----------------", iter);
#ifndef ERASE_CENTER
		for (list<st_center*>::iterator it = center_list.begin(); it != center_list.end(); it++) {
			scan_center(*it, elem_list.begin(), elem_list.end());
		}
#else
		for (list<st_center*>::iterator it = center_list.begin(); it != center_list.end(); ) {
			if ((*it)->one->similar >= SPECIES_SOFT) {
				center_list.erase(it++);
			}else{
				scan_center(*it, elem_list.begin(), elem_list.end());
				it++;
			}
		}
#endif
		for (list<st_elem*>::iterator it = elem_list.begin(); it != elem_list.end(); it++) {
			st_elem* p = *it;
			if (p->ctr && p->ctr->one == p) { 
				//all has been done in scan_center.
				//logger->info("cluster: all has been done in scan_center.");
			}
			else if(p->similar >= UNION_THRE) {
				p->ctr->one->union_elem(p);
				//logger->info("cluster: union");
			}
			else if (p->similar >= SPECIES_SOFT) {
				p->ctr->grp.push_back(p);
				//logger->info("cluster: cluster");

			}
			else{
				//logger->info("cluster: new, nearest = %.4f", p->similar);
				center_list.push_back(new st_center(p));
				scan_center(center_list.back(), it, elem_list.end());
			}
			//p->ctr = NULL;
		}

		logger->info("-------------------ITER=[%d] scaned, timecost=[%.2lf]s------", iter, COUNT_SECONDS(start));
		if (iter) {
			start = clock();
			int active_size = center_list.size();
			for (list<st_center*>::iterator it = center_list.begin(); it != center_list.end();)
				if ((*it)->recenter(align, SPECIES_SOFT)) {
					fix_list.push_back(*it);
					center_list.erase(it++);

					/*
					   logger->info("debug: move");
					   st_center *p = fix_list.back();
					   for (list<st_elem*>::iterator jt = p->grp.begin(); jt != p->grp.end(); jt++) {
					   printf("%d ", (*jt)->n);
					   for (list<int>::iterator kt = (*jt)->union_list.begin(); kt != (*jt)->union_list.end(); kt++)
					   printf("%d ", (*kt));
					   }
					   printf("\n");
					   */

				}else{
					it++;

				}
			logger->info("move to fixed=[%d], cluster size=[%d], recenter timecost=[%.2lf]s", active_size - center_list.size(), center_list.size() + fix_list.size(), COUNT_SECONDS(start));
		}

		start = clock();
		for (list<st_elem*>::iterator it = elem_list.begin(); it != elem_list.end();) {
			switch((*it)->status) {
				case st_elem::UNIONED: delete (*it); elem_list.erase(it++); break;
				case st_elem::FIXED  : elem_list.erase(it++); break;
				case st_elem::ACTIVE : (*it)->ctr = NULL; (*it)->similar = 0; it++; break;
				default: logger->info("ERROR: elem unknown status");
			}
		}

/*
#ifdef __DEBUG__
		int s1 = 0, s2 =0;
		for (list<st_center*>::iterator it = center_list.begin(); it != center_list.end(); it++)
			if ((*it)->grp.size() != 1)
				printf("center=%d, grpsize=%ld\n", (*it)->one->n, (*it)->grp.size());
			else
				s1++;

		for (list<st_center*>::iterator it = fix_list.begin(); it != fix_list.end(); it++)
			if ((*it)->grp.size() != 1)
				printf("solid_center=%d, grpsize=%ld\n", (*it)->one->n, (*it)->grp.size());
			else
				s2++;

		printf("ONE ELEM cluster: %d\n", s1);
		printf("ONE ELEM solid_cluster: %d\n", s2);
		printf("TOTAL cluster num= %ld\n", fix_list.size() + center_list.size());
		printf("-------------------END ITER=%d-----------------\n\n", iter);
#endif
*/
	}
}

void cls_cluster::print_result() {
	for (list<st_center*>::iterator it = center_list.begin(); it != center_list.end(); it++) {
		st_center* p = *it;
		for (vector<st_elem*>::iterator jt = p->grp.begin(); jt != p->grp.end(); jt++) {
			printf("%d ", (*jt)->n);
			for (list<int>::iterator kt = (*jt)->union_list.begin(); kt != (*jt)->union_list.end(); kt++)
				printf("%d ", (*kt));
		}
		printf("\n");
	}
		
	for (list<st_center*>::iterator it = fix_list.begin(); it != fix_list.end(); it++) {
		st_center* p = *it;
		for (vector<st_elem*>::iterator jt = p->grp.begin(); jt != p->grp.end(); jt++) {
			printf("%d ", (*jt)->n);
			for (list<int>::iterator kt = (*jt)->union_list.begin(); kt != (*jt)->union_list.end(); kt++)
				printf("%d ", (*kt));
		}
		printf("\n");
	}
}


void cls_cluster::export_result(vector<vector<long> >& result, vector<int>& weight) 
{
	for (list<st_center*>::iterator it = center_list.begin(); it != center_list.end(); it++) {
		st_center* p = *it;
		result.push_back(vector<long>());
		weight.push_back(0);
		for (vector<st_elem*>::iterator jt = p->grp.begin(); jt != p->grp.end(); jt++) {
			result.back().push_back((*jt)->n);
			weight.back() += (*jt)->weight;
			for (list<int>::iterator kt = (*jt)->union_list.begin(); kt != (*jt)->union_list.end(); kt++)
				result.back().push_back(*kt);
		}
	}
		
	for (list<st_center*>::iterator it = fix_list.begin(); it != fix_list.end(); it++) {
		st_center* p = *it;
		result.push_back(vector<long>());
		weight.push_back(0);
		for (vector<st_elem*>::iterator jt = p->grp.begin(); jt != p->grp.end(); jt++) {
			result.back().push_back((*jt)->n);
			weight.back() += (*jt)->weight;
			for (list<int>::iterator kt = (*jt)->union_list.begin(); kt != (*jt)->union_list.end(); kt++)
				result.back().push_back(*kt);
		}
	}
}
