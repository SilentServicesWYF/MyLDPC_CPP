#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;
void readdata(int row,int col, int **matrix, char const *filename)
{
    ifstream fin;
    fin.open(filename);
    for (int i = 0; i < row; i++)
    {    
        for (int j = 0; j < col; j++)
        {
            fin >> matrix[i][j];
        }
    }
    std::cout<<filename<<"矩阵导入完成"<<endl;
}

void readvector(int row, int *vector, char const *filename)
{
    ifstream fin;
    fin.open(filename);
    for (int i = 0; i < row; i++)
    {
        fin >> vector[i];
    }
    std::cout<<filename<<"向量导入完成"<<endl;
}

int *nodediff(int *mset, int m, int n2m_num)
{
    int *msubset = new int[n2m_num-1];
    int msubset_index = 0;
    for (int k = 0; k < n2m_num; k ++)
    {
        if (mset[k] != m)
        {
            msubset[msubset_index] = mset[k];
            msubset_index = msubset_index + 1;
        }
    }
    return msubset;
}

float *floatslice(float *arry, int start_index, int end_index)
/*用于float类型的切片*/
{
    float *a = new float[end_index-start_index + 1];
    int count = 0;
    for (int i = start_index; i <= end_index;i++)
    {
        a[count] = arry[i];
        count = count + 1;
    }
    return a;
}

int *intslice(int *arry, int start_index, int end_index)
/*用于int类型的切片*/
{
    int *a = new int[end_index-start_index];
    int count = 0;
    for (int i = start_index; i <= end_index;i++)
    {
        a[count] = arry[i];
        count = count + 1;
    }
    return a;
}

int findnode(int **arry, int node1, int node2,int len)
{
    int *temp;
    temp = arry[node1];
    for (int k = 0; k < len; k++)
    {
        if (temp[k] == node2)
        {
            return k;
            break;
        }
    }
}