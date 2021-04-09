#ifndef AUXIL_H
#define AUXIL_H

void readdata(int row,int col, int **matrix, char const *filename);
float *floatslice(float *arry, int start_index, int end_index);
int *intslice(int *arry, int start_index, int end_index);
int findnode(int **arry, int node1, int node2,int len);
void readvector(int row, int *vector, char const *filename);
int *nodediff(int *mset, int m, int n2m_num);

#endif