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
#include <stdarg.h>
#include <ctime>
#include <cstdio>

#include "cls_log.h"

using namespace std;

cls_log *cls_log::instance(NULL);

cls_log::cls_log() {}

cls_log::~cls_log() {}

cls_log* cls_log::get_instance()
{
	if (instance == NULL)
		instance = new cls_log();
	return instance;
}

void cls_log::set_mode(int _mode)
{
	mode = _mode;
}


void cls_log::info(const char* format, ...)
{
	time_t lt = time(NULL);
	struct tm *ptr = localtime(&lt);
	char buf[80];
	strftime(buf, 80, "%X", ptr);

	int MPI_RANK;
	MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);

	if (MPI_RANK == 0) {
		fprintf(stdout, "[%s][INFO][Master] ", buf);
	}else if (mode == 1) {
		fprintf(stdout, "[%s][INFO][Rank %d] ", buf, MPI_RANK);
	}else {
		return ;
	}
	
	/*
	fprintf(stdout, "[%s]: ", buf);
	*/
	va_list args;
	va_start(args, format);
	vfprintf(stdout, format, args);
	va_end(args);
	fprintf(stdout, "\n");
}

void cls_log::fatal(const char* format, ...)
{
	time_t lt = time(NULL);
	struct tm *ptr = localtime(&lt);
	char buf[80];
	strftime(buf, 80, "%X", ptr);

	int MPI_RANK;
	MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);

	if (MPI_RANK == 0) {
		fprintf(stdout, "[%s][FATAL][Master] ", buf);
	}else if (mode == 1) {
		fprintf(stdout, "[%s][FATAL][Rank %d] ", buf, MPI_RANK);
	}else {
		return ;
	}
	
	/*
	fprintf(stdout, "[%s]: ", buf);
	*/
	va_list args;
	va_start(args, format);
	vfprintf(stdout, format, args);
	va_end(args);
	fprintf(stdout, "\n");
}
