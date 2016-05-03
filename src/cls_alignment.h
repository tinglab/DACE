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
 * 
 *   *************** NOTE ****************
 *      We are using the source code from project CD-HIT (written by Dr. Weizhong Li, http://cd-hit.org/)
 *      under GPLv2 license in the alignment function, in order to make the clustering results consistent
 *      with the state-of-the-art programs. 
 *
 */


#ifndef __cls__alignment__H__
#define __cls__alignment__H__

#include <string>
#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

class cls_alignment
{
public:

	enum { DP_BACK_NONE=0, DP_BACK_LEFT_TOP=1, DP_BACK_LEFT=2, DP_BACK_TOP=3 };

	static const unsigned int HASH_SIZE = 1048576 * 8 - 1;   //2^23-1

	cls_alignment(const std::vector<std::string>& raw, bool _using_fullhash = false);
	~cls_alignment();

	double get_similarity(long nx, long ny, double threshold);
	double get_similarity(int nx, int ny, double threshold);

private:
	
	int diag_test_aapn_est(char iseq1[], char iseq2[], int len1, int len2, int &best_sum, int band_width, int &band_left, int &band_center, int &band_right, double cluster_thd, bool recalc_aap);

	int local_band_align(char iseq1[], char iseq2[], int len1, int len2, int &best_score, int &iden_no, int &alnln, float &dist, int *alninfo, int band_left, int band_center, int band_right);

	unsigned char char_conv(char ch) ;


private:
	
	int taap[256], aap_begin[256];
	int *aap_list, *diag_score, *diag_score2;
	int **score_mat, **back_mat;

	int *dist_hash_n;
	double *dist_hash_v;
	boost::unordered_map<std::pair<size_t, size_t>, double > fullhash;

	char *seq1, *seq2;

	bool using_fullhash;
	const std::vector<std::string>& raw;
	long cnt1, cnt2, lastnx;
};

#endif
