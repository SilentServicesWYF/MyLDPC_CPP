#include "gfcalu.h"
#include <stdio.h>

int gfaddtable[4][4] = {{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
int gfsubtable[4][4] = {{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
int gfmultable[4][4] = {{0,0,0,0},{0,1,2,3},{0,2,3,1},{0,3,1,2}};
int gfdivtable[4][4] = {{0,0,0,0},{0,1,3,2},{0,2,1,3},{0,3,2,1}};

int gfadd(int a, int b)
{
    return gfaddtable[a][b];
}
int gfsub(int a, int b)
{
    return gfsubtable[a][b];
}
int gfmul(int a, int b)
{
    return gfmultable[a][b];
}
int gfdiv(int a, int b)
{
    return gfdivtable[a][b];
}
int **gfmatrixmul(int **mat1, int **mat2, int row1, int col1)
/*验证H*c的函数*/
{
    int **mat = new int *[row1];
    for (int m = 0; m < row1; m++)
    {
        mat[m] = new int[1];
        mat[m][0] = 0;
        for (int n = 0;n < col1; n++)
        {
            mat[m][0] = gfadd(mat[m][0],gfmul(mat1[m][n],mat2[n][0]));
        }
    }
    return mat;
}