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

#include <cstring>
#include <string>
#include <iostream>

#include "cls_alignment.h"
#include "cls_common.h"

using namespace std;
#include <cstdlib>
#include <cstdio>

#define MAX_SEQ_LEN 2048

unsigned char cls_alignment::char_conv(char ch) 
{
	switch (ch) {
		case 'A': case 'a': return 0;
		case 'T': case 't': return 1;
		case 'C': case 'c': return 2;
		case 'G': case 'g': return 3;
		default:
			return 4;
	}
}


int cls_alignment::diag_test_aapn_est(char iseq1[], char iseq2[], int len1, int len2, int &best_sum, int band_width, int &band_left, int &band_center, int &band_right, double cluster_thd, bool recalc_aap)
{
	int i, i1, j, k;
	int *pp, *pp2;
	int nall = len1+len2-1;
	int NAA1 = 4;
	int NAA2 = NAA1 * NAA1;
	int NAA3 = NAA2 * NAA1;

/*	static int taap[256], aap_begin[256];
	static int aap_list[MAX_SEQ_LEN * 2];
	static int diag_score [MAX_SEQ_LEN * 2];
	static int diag_score2[MAX_SEQ_LEN * 2];
*/
	//if (nall > MAX_DIAG) bomb_error("in diag_test_aapn_est, MAX_DIAG reached");
	pp = & diag_score[0];
	pp2 = & diag_score2[0];
	
	memset(diag_score, 0, nall * sizeof(int));
	memset(diag_score2, 0, nall * sizeof(int));

if (recalc_aap) {
	memset(taap, 0, 256 * sizeof(int));
	for (i = 0; i < len1 - 3; i++) {
		unsigned char c0 = (iseq1[i+0]);
		unsigned char c1 = (iseq1[i+1]);
		unsigned char c2 = (iseq1[i+2]);
		unsigned char c3 = (iseq1[i+3]);
		if ((c0>=4) || (c1>=4) || (c2>=4) || (c3>=4)) continue; //skip N

		int c22 = c0*NAA3+ c1*NAA2 + c2*NAA1 + c3;
		taap[c22]++;
	}
	aap_begin[0] = taap[0];
	for (i = 1; i < 256; i++)
		aap_begin[i] = aap_begin[i-1] + taap[i];
	for (i = 0; i < len1 - 3; i++) {
		unsigned char c0 = (iseq1[i+0]);
		unsigned char c1 = (iseq1[i+1]);
		unsigned char c2 = (iseq1[i+2]);
		unsigned char c3 = (iseq1[i+3]);
		if ((c0>=4) || (c1>=4) || (c2>=4) || (c3>=4)) continue; //skip N

		int c22 = c0*NAA3+ c1*NAA2 + c2*NAA1 + c3;
		aap_list[--aap_begin[c22]] = i;
	}
}

	int *bip;
	int c22, cpx;
	int len22 = len2-3;
	i1 = len1-1;
	for (i=0; i<len22; i++,i1++,iseq2++) {
		unsigned char c0 = (iseq2[0]);
		unsigned char c1 = (iseq2[1]);
		unsigned char c2 = (iseq2[2]);
		unsigned char c3 = (iseq2[3]);
		if ((c0>=4) || (c1>=4) || (c2>=4) || (c3>=4)) continue; //skip N

		c22 = c0*NAA3+ c1*NAA2 + c2*NAA1 + c3;
		if ( (j=taap[c22]) == 0) continue;
		cpx = 1 + (c0 != c1) + (c1 != c2) + (c2 != c3);
		bip = & aap_list[ aap_begin[c22] ];	 //	bi = aap_begin[c22];
		for (; j; j--, bip++) { 
			diag_score[i1 - *bip]++;
			diag_score2[i1 - *bip] += cpx;
		}
	}

	//find the best band range
	int required_aa1 = int(cluster_thd * (double) len2 - 3);
	int band_b = max(required_aa1-1, 0);  // on dec 21 2001
	int band_e = nall - required_aa1;

	int band_m = min( band_b+band_width-1, band_e );
	int best_score=0;
	int best_score2=0;
	int max_diag = 0;
	int max_diag2 = 0;
	int imax_diag = 0;
	for (i=band_b; i<=band_m; i++){
		best_score += diag_score[i];
		best_score2 += diag_score2[i];
		if( diag_score2[i] > max_diag2 ){
			max_diag2 = diag_score2[i];
			max_diag = diag_score[i];
			imax_diag = i;
		}
	}
	int from=band_b;
	int end =band_m;
	int score = best_score;  
	int score2 = best_score2;  

	for (k=from, j=band_m+1; j<band_e; j++, k++) {
		score -= diag_score[k]; 
		score += diag_score[j]; 
		score2 -= diag_score2[k]; 
		score2 += diag_score2[j]; 
		if ( score2 > best_score2 ) {
			from = k + 1;
			end  = j;
			best_score = score;
			best_score2 = score2;
			if( diag_score2[j] > max_diag2 ){
				max_diag2 = diag_score2[j];
				max_diag = diag_score[j];
				imax_diag = j;
			}
		}
	}
#if 0
	printf( "%i\n", required_aa1 );
	printf( "max=%3i  imax=%3i; band:  %3i  %3i  %i\n", max_diag, imax_diag, band_b, band_e, band_m );
	printf( "best: %i\n", best_score );
	printf( "from: %i, end: %i,  best: %i\n", from, end, best_score );
#endif
	int mlen = imax_diag;
	if( imax_diag > len1 ) mlen = nall - imax_diag;
	int emax = int((1.0 - cluster_thd) * mlen) + 1;
	for (j=from; j<imax_diag; j++) { // if aap pairs fail to open gap
		if ( (imax_diag - j) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; from++;
		} else break;
	}
	for (j=end; j>imax_diag; j--) { // if aap pairs fail to open gap
		if ( (j - imax_diag) > emax || diag_score[j] < 1 ) {
			best_score -= diag_score[j]; end--;
		} else break;
	}
	//printf( "best: %i\n", best_score );

	band_left = from-len1+1; 
	band_right= end-len1+1;
	band_center = imax_diag - len1 + 1;
	best_sum = best_score;
#if 0
	printf( "max=%3i  imax=%3i; band:  %3i  %3i  %i\n", mmax, immax, band_b, band_e, band_m );
	printf( "%3i:  best: %i,  %i  %i  %i\n", required_aa1, best_score, band_left, band_right, band_width );
#endif
	return 0;
}


int cls_alignment::local_band_align(char iseq1[], char iseq2[], int len1, int len2, int &best_score, int &iden_no, int &alnln, float &dist, int *alninfo, int band_left, int band_center, int band_right)
{

/*
	printf( "int len1=%d, int len2=%d, \
		int &best_score=%d, int &iden_no=%d, int &alnln=%d, float &dist=%f, int *alninfo=%d,\
		int band_left=%d, int band_center=%d, int band_right=%d)",
		len1, len2, 
		best_score, iden_no, alnln, dist, *alninfo,
		band_left, band_center, band_right);
*/

	int mat_gap = -3932160, mat_ext_gap = -655360;
	int mat_matrix[4][4];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mat_matrix[i][j] = 1310720 * (i == j ? 1 : -1);

	int i, j, k, j1;
	int jj, kk;
	int iden_no1;
	int64_t best_score1;
	iden_no = 0;

	if ( (band_right >= len2 ) ||
			(band_left  <= -len1) ||
			(band_left  > band_right) ) {
		//printf("fuck error\n");
		return -1;
	}

	// allocate mem for score_mat[len1][len2] etc
	int band_width = band_right - band_left + 1;
	int band_width1 = band_width + 1;

	//printf( "%i  %i\n", band_right, band_left );

/*	static int score_mat[MAX_SEQ_LEN][MAX_SEQ_LEN];
	static int back_mat[MAX_SEQ_LEN][MAX_SEQ_LEN];
*/

	best_score = 0;
	/*
	   seq2 len2 = 17			seq2 len2 = 17	  seq2 len2 = 17
	   01234567890123456	   01234567890123456	01234567890123456
	   0	 xxxxxxxxxxxxxxxxx \\\\\\XXXxxxxxxxxxxxxxx	xXXXXXXXxxxxxxxxx
	   1\\\\\Xxxxxxxxxxxxxxxxx  \\\\\Xxx\xxxxxxxxxxxxx	xx\xxxxx\xxxxxxxx
	   2 \\\\X\xxxxxxxxxxxxxxx   \\\\Xxxx\xxxxxxxxxxxx	xxx\xxxxx\xxxxxxx
	   seq1 3  \\\Xx\xxxxxxxxxxxxxx	\\\Xxxxx\xxxxxxxxxxx	xxxx\xxxxx\xxxxxx
	   len1 4   \\Xxx\xxxxxxxxxxxxx	 \\Xxxxxx\xxxxxxxxxx	xxxxx\xxxxx\xxxxx
	   = 11 5	\Xxxx\xxxxxxxxxxxx	  \Xxxxxxx\xxxxxxxxx	xxxxxx\xxxxx\xxxx
	   6	 Xxxxx\xxxxxxxxxxx	   Xxxxxxxx\xxxxxxxx	xxxxxxx\xxxxx\xxx
	   7	 x\xxxx\xxxxxxxxxx	   x\xxxxxxx\xxxxxxx	xxxxxxxx\xxxxx\xx
	   8	 xx\xxxx\xxxxxxxxx	   xx\xxxxxxx\xxxxxx	xxxxxxxxx\xxxxx\x
	   9	 xxx\xxxx\xxxxxxxx	   xxx\xxxxxxx\xxxxx	xxxxxxxxxx\xxxxx\
	   0	 xxxx\xxxx\xxxxxxx	   xxxx\xxxxxxx\xxxx	xxxxxxxxxxx\xxxxx
	   band_left < 0		   band_left < 0		band_left >=0
	   band_right < 0		  band_right >=0	   band_right >=0
	   init score_mat, and iden_mat (place with upper 'X')
	 */

	if (band_left < 0) {  //set score to left border of the matrix within band
		int tband = (band_right < 0) ? band_right : 0;
		//for (k=band_left; k<tband; k++)
		for (k=band_left; k<=tband; k++) { // fixed on 2006 11 14
			i = -k;
			j1 = k-band_left;
			// penalty for leading gap opening = penalty for gap extension
			score_mat[i][j1] =  mat_ext_gap * i;
			back_mat[i][j1] = DP_BACK_TOP;
		}
		back_mat[-tband][tband-band_left] = DP_BACK_NONE;
	}

	if (band_right >=0) { //set score to top border of the matrix within band
		int tband = (band_left > 0) ? band_left : 0;
		for (j=tband; j<=band_right; j++) {
			j1 = j-band_left;
			score_mat[0][j1] = mat_ext_gap * j;
			back_mat[0][j1] = DP_BACK_LEFT;
		}
		back_mat[0][tband-band_left] = DP_BACK_NONE;
	}

	int gap_open[2] = { mat_gap, mat_ext_gap };
	int max_diag = band_center - band_left;
	int extra_score[4] = { 4, 3, 2, 1 };
	for (i=1; i<=len1; i++) {
		int J0 = 1 - band_left - i;
		int J1 = len2 - band_left - i;
		if( J0 < 0 ) J0 = 0;
		if( J1 >= band_width ) J1 = band_width;
		for (j1=J0; j1<=J1; j1++){
			j = j1+i+band_left;

			int ci = iseq1[i-1];
			int cj = iseq2[j-1];
			int sij = mat_matrix[ci][cj];
			//int iden_ij = (ci == cj);
			int s1, k0, back;

			/* extra score according to the distance to the best diagonal */
			int extra = extra_score[ abs(j1 - max_diag) & 3 ]; // max distance 3
			sij += extra * (sij>0);

			back = DP_BACK_LEFT_TOP;
			best_score1 = score_mat[i-1][j1] + sij;
			int gap0 = gap_open[ (i == len1) | (j == len2) ];
			int gap = 0;
			int64_t score;

			if( j1 > 0 ){
				gap = gap0;
				if( back_mat[i][j1-1] == DP_BACK_LEFT ) gap = mat_ext_gap;
				if( (score = score_mat[i][j1-1] + gap) > best_score1 ){
					back = DP_BACK_LEFT;
					best_score1 = score;
				}
			}
			if(j1+1<band_width){
				gap = gap0;
				if( back_mat[i-1][j1+1] == DP_BACK_TOP ) gap = mat_ext_gap;
				if( (score = score_mat[i-1][j1+1] + gap) > best_score1 ){
					back = DP_BACK_TOP;
					best_score1 = score;
				}
			}
			score_mat[i][j1] = best_score1;
			back_mat[i][j1]  = back;
			//printf( "%2i(%2i) ", best_score1, iden_no1 );

		}
		//printf( "\n" );
	}
	i = j = 0;
	if( len2 - band_left < len1 ){
		i = len2 - band_left;
		j = len2;
	}else if( len1 + band_right < len2 ){
		i = len1;
		j = len1 + band_right;
	}else{
		i = len1;
		j = len2;
	}
	j1 = j - i - band_left;
	best_score = score_mat[i][j1];
	best_score1 = score_mat[i][j1];

#if 1
	const char *letters = "acgtnx";
	const char *letters2 = "ACGTNX";
#else
	const char *letters = "arndcqeghilkmfpstwyvbzx";
	const char *letters2 = "ARNDCQEGHILKMFPSTWYVBZX";
#endif
	int back = back_mat[i][j1];
	int last = back;
	int count = 0, count2 = 0, count3 = 0;
	int match, begin1, begin2, end1, end2;
	int gbegin1=0, gbegin2=0, gend1=0, gend2=0;
	int64_t score, smin = best_score1, smax = best_score1 - 1;
	int posmin, posmax, pos = 0;
	int bl, dlen = 0, dcount = 0;
	posmin = posmax = 0;
	begin1 = begin2 = end1 = end2 = 0;

#ifdef PRINT
#define PRINT
	printf( "%i %i\n", best_score, score_mat[i][j1] );
	printf( "%i %i %i\n", band_left, band_center, band_right );
	printf( "%i %i %i %i\n", i, j, j1, len2 );
#endif
	int masked = 0;
	int indels = 0;
	int max_indels = 0;

	bool global_identity = true;
	while( back != DP_BACK_NONE ){
		switch( back ){
		case DP_BACK_TOP  :
#ifdef PRINT
			printf( "%5i: %c %c %9i\n", pos, letters[ iseq1[i-1] ], '|', score_mat[i][j1] );
#endif
			bl = (last != back) & (j != 1) & (j != len2);
			dlen += bl;
			dcount += bl;
			score = score_mat[i][j1];
			if( score < smin ){
				count2 = 0;
				smin = score;
				posmin = pos - 1;
				begin1 = i;
				begin2 = j;
			}
			i -= 1;
			j1 += 1;
			break;
		case DP_BACK_LEFT :
#ifdef PRINT
			printf( "%5i: %c %c %9i\n", pos, '|', letters[ iseq2[j-1] ], score_mat[i][j1] );
#endif
			bl = (last != back) & (i != 1) & (i != len1);
			dlen += bl;
			dcount += bl;
			score = score_mat[i][j1];
			if( score < smin ){
				count2 = 0;
				smin = score;
				posmin = pos - 1;
				begin1 = i;
				begin2 = j;
			}
			j1 -= 1;
			j -= 1;
			break;
		case DP_BACK_LEFT_TOP :
			if( alninfo && global_identity ){
				if( i == 1 || j == 1 ){
					gbegin1 = i-1;
					gbegin2 = j-1;
				}else if( i == len1 || j == len2 ){
					gend1 = i-1;
					gend2 = j-1;
				}
			}
			score = score_mat[i][j1];
			i -= 1;
			j -= 1;
			match = iseq1[i] == iseq2[j];
			if( score > smax ){
				count = 0;
				smax = score;
				posmax = pos;
				end1 = i;
				end2 = j;
			}
			dlen += 1;
			dcount += ! match;
			count += match;
			count2 += match;
			count3 += match;

			if( score < smin ){
				int mm = match == 0;
				count2 = 0;
				smin = score;
				posmin = pos - mm;
				begin1 = i + mm;
				begin2 = j + mm;
			}
			break;
		default : printf( "%i\n", back ); break;
		}
		pos += 1;
		last = back;
		back = back_mat[i][j1];
	}
	iden_no = global_identity ? count3 : count - count2;
	alnln = posmin - posmax + 1 - masked;
	dist = dcount/(float)dlen;
	//dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
	int umtail1 = len1 - 1 - end1;
	int umtail2 = len2 - 1 - end2;
	int umhead = begin1 < begin2 ? begin1 : begin2;
	int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
	int umlen = umhead + umtail;
	
	if( alninfo ){
		alninfo[0] = begin1;
		alninfo[1] = end1;
		alninfo[2] = begin2;
		alninfo[3] = end2;
		alninfo[4] = masked;
		if( global_identity ){
			alninfo[0] = gbegin1;
			alninfo[1] = gend1;
			alninfo[2] = gbegin2;
			alninfo[3] = gend2;
		}
	}
#ifdef PRINT
	printf( "%6i %6i:  %4i %4i %4i %4i\n", alnln, iden_no, begin1, end1, begin2, end2 );
	printf( "%6i %6i:  %4i %4i\n", posmin, posmax, posmin - posmax, count - count2 );
	printf( "smin = %9i, smax = %9i\n", smin, smax );
	printf( "dlen = %5i, dcount = %5i, dist = %.5f\n", dlen, dcount, dcount/(float)dlen );
#endif
	
	return 0;
} 


double cls_alignment::get_similarity(long nx, long ny, double threshold)
{
	if (nx == ny)
		return 1.0;

	long n1 = max(nx, ny), n2 = min(nx, ny);

	pair<size_t, size_t> pr = make_pair(n1, n2);

	if (fullhash.find(pr) != fullhash.end()) {
		cnt2 += 1;
		return fullhash[pr];
	}

	double result;
	cnt1 += 1;

	if (raw[nx].size() < raw[ny].size()) swap(nx, ny);
	
	for (int j = 0; j < raw[nx].size(); j++) seq1[j] = char_conv(raw[nx][j]);
	for (int j = 0; j < raw[ny].size(); j++) seq2[j] = char_conv(raw[ny][j]);
	
	int band_width = 20;
	int best_sum, band_left, band_center, band_right;
	
	int OK = diag_test_aapn_est(seq1, seq2, raw[nx].size(), raw[ny].size(), best_sum, band_width, band_left,  band_center, band_right, threshold, lastnx != nx);
	lastnx = nx;

	if (OK == 0) {
		//printf("best_sum %d, band_left %d, band_center %d, band_right %d \n",  best_sum, band_left, band_center, band_right);

		int best_score, iden_no, alnln;
		float dist;
		OK = local_band_align(seq1, seq2, raw[nx].size(), raw[ny].size(), best_score, iden_no, alnln, dist, NULL,  band_left,  band_center,  band_right);
		
		if (OK == -1) {
			dist = 0.5;
		}
		result = double(1 - dist);

	}else{
		result = 0.5;
	}
	
	fullhash[pr] = result;

	return result;
}
double cls_alignment::get_similarity(int nx, int ny, double threshold)
{
	if (nx == ny)
		return 1.0;

	int n1 = max(nx, ny),
	    n2 = min(nx, ny),
	    hash_pos = n1 * (n2+1) % HASH_SIZE;

	if (dist_hash_n[hash_pos] == n1) {
		cnt2 += 1;
		return dist_hash_v[hash_pos];
	}

	double result;
	cnt1 += 1;

	if (raw[nx].size() < raw[ny].size()) swap(nx, ny);
	
	for (int j = 0; j < raw[nx].size(); j++) seq1[j] = char_conv(raw[nx][j]);
	for (int j = 0; j < raw[ny].size(); j++) seq2[j] = char_conv(raw[ny][j]);
	
	int band_width = 20;
	int best_sum, band_left, band_center, band_right;
	
	int OK = diag_test_aapn_est(seq1, seq2, raw[nx].size(), raw[ny].size(), best_sum, band_width, band_left,  band_center, band_right, threshold, lastnx != nx);
	lastnx = nx;

	if (OK == 0) {
		//printf("best_sum %d, band_left %d, band_center %d, band_right %d \n",  best_sum, band_left, band_center, band_right);

		int best_score, iden_no, alnln;
		float dist;
		OK = local_band_align(seq1, seq2, raw[nx].size(), raw[ny].size(), best_score, iden_no, alnln, dist, NULL,  band_left,  band_center,  band_right);
		
		if (OK == -1) {
			dist = 0.5;
		}
		result = double(1 - dist);

	}else{
		result = 0.5;
	}
	
	dist_hash_n[hash_pos] = n1;
	dist_hash_v[hash_pos] = result;

	return result;
}


cls_alignment::cls_alignment(const vector<string>& raw, bool _using_fullhash)
	: raw(raw), cnt1(0), cnt2(0), lastnx(-1), using_fullhash(_using_fullhash)
{

	if (!_using_fullhash) {
		dist_hash_n = new int [HASH_SIZE];
		dist_hash_v = new double [HASH_SIZE];
		memset(dist_hash_n, 0xFFFF, sizeof(int) * HASH_SIZE);
	}

	seq1 = new char[MAX_SEQ_LEN];
	seq2 = new char[MAX_SEQ_LEN];


	aap_list = new int [MAX_SEQ_LEN * 2];
	diag_score = new int [MAX_SEQ_LEN * 2];
	diag_score2 = new int[MAX_SEQ_LEN * 2];
	
	score_mat = create_array<int>(MAX_SEQ_LEN, MAX_SEQ_LEN);
	back_mat = create_array<int>(MAX_SEQ_LEN, MAX_SEQ_LEN);
}

cls_alignment::~cls_alignment()
{
	//printf("cls_alignment.count = %ld, cls_alignemnt.hash_count = %ld, fullhash = %ld\n", cnt1, cnt2, fullhash.size());
	if (!using_fullhash) {
		delete [] dist_hash_n;
		delete [] dist_hash_v;
	}

	delete [] seq1;
	delete [] seq2;

	delete [] aap_list;
	delete [] diag_score;
	delete [] diag_score2;

	delete_array(score_mat, MAX_SEQ_LEN, MAX_SEQ_LEN);
	delete_array(back_mat, MAX_SEQ_LEN, MAX_SEQ_LEN);
}


