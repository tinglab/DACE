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

#ifndef __cls__log__H__
#define __cls__log__H__

class cls_log
{

private:
	cls_log();
	~cls_log();

public:
	static cls_log* get_instance();

	void set_mode(int mode);

	void info(const char* format, ...);
	void fatal(const char* format, ...);

private:
	static cls_log *instance;

	int mode;
};

#endif

