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

#ifndef __cls__common__H__
#define __cls__common__H__

#include <cstdlib>

//#define __DEBUG__LOG__

#define COUNT_SECONDS(start) (double(clock() - start) / CLOCKS_PER_SEC)

template <class T> const T& max2 (const T& a, const T& b) {
  return (a < b) ? b : a; 
}

template <class T> const T& min2 (const T& a, const T& b) {
  return (a > b) ? b : a; 
}

template <class T> const T& max3 (const T& a, const T& b, const T& c) {
  return (a > b) ? max2(a, c) : max2(b, c);
}

template <class T> const T& min3 (const T& a, const T& b, const T& c) {
  return (a < b) ? min2(a, c) : min2(b, c);
}


template <class T> T** create_array(int row, int col)
{
	int size = sizeof(T);
	int point_size = sizeof(T*);
	T **arr = (T **) malloc(point_size * row + size * row * col);
	if (arr != NULL)
	{	
		T *head = (T*)((long)arr + point_size * row);
		for (int i = 0; i < row; ++i)
		{
			arr[i] =  (T*)((long)head + i * col * size);
			for (int j = 0; j < col; ++j)
				new (&arr[i][j]) T;
		}
	}
	//debug("create_array: size=[%d]x[%d]", row, col);
	return (T**)arr;
}

template <class T> void delete_array(T **arr, int row, int col)
{
	for (int i = 0; i < row; ++i)
		for (int j = 0; j < col; ++j)
			arr[i][j].~T();
	if (arr != NULL)
		free((void**)arr);
}

#endif
