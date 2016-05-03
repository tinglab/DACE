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


#ifndef __cls__kmer__H__
#define __cls__kmer__H__

class cls_kmer
{

/*
  Use 2 bits for each base pair:
	array[size]: ....
	 .
	 .
	 .
	array[1]:  s[63] s[62] ... s[32]
	array[0]:  s[31] s[30] ... s[0]
*/

public:
	cls_kmer(int length, const char* s);
	cls_kmer(const cls_kmer &b);
	~cls_kmer();

	void push_back(char c);
	void print() const;
	bool operator< (const cls_kmer &b) const;
	long operator% (const long m) const;
	bool operator== (const cls_kmer &b) const;
	
	struct hash_value {
		unsigned long operator() (const cls_kmer & a) const;
	};
	
private:
	int kmerlen, size;
	unsigned long *array;
};

#endif

