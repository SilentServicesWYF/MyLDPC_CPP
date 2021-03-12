#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;
int gfaddtable[4][4] = {{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
int gfsubtable[4][4] = {{0,1,2,3},{1,0,3,2},{2,3,0,1},{3,2,1,0}};
int gfmultable[4][4] = {{0,0,0,0},{0,1,2,3},{0,2,3,1},{0,3,1,2}};
int gfdivtable[4][4] = {{0,0,0,0},{0,1,3,2},{0,2,1,3},{0,3,2,1}};

int gfadd(int a, int b);
int gfsub(int a, int b);
int gfmul(int a, int b);
int gfdiv(int a, int b);
int **gfmatrixmul(int **mat1, int **mat2, int row1, int col1);
void readdata(int row,int col, int **matrix, char *filename);
float LLRV(float *subconstell, int *pskdict, int gf);
float *floatslice(float *arry, int start_index, int end_index);
int *intslice(int *arry, int start_index, int end_index);

int main()
{
    //参数设定
    int row1 = 186;
    int col1 = 372;
    int row2 = 372;
    int col2 = 1;
    int maxweight1 = 3;
    int maxweight2 = 5;
    //初始化内存
    int ** H = new int *[row1];
    int ** c = new int *[row2];
    int ** flag = new int *[row1]; //flag在每次循环中判断跳出条件
    int ** m2n = new int *[row1]; //与m连接的n的连接表
    int ** n2m = new int *[row2]; //与n连接的m的连接表
    float constell[2*row2] = {0};
    float *subconstell;
    for (int i = 0; i<row1; i++)
    {        
        flag[i] = new int[1];
        H[i] = new int[col1];
        m2n[i] = new int[maxweight2];
    }
    for (int i = 0; i<row2; i++)
    {
        n2m[i] = new int[maxweight1];
        c[i] = new int[col2];
    }
    //读取文件
    readdata(row1,col1,H,"H.txt");
    readdata(row2,col2,c,"c.txt");
    readdata(row1,maxweight2,m2n,"m2n.txt");
    readdata(row2,maxweight1,n2m,"n2m.txt");
    ifstream fin;
    fin.open("constell.txt");
    for (int i = 0; i < 2*row2; i++)
    {
        fin >> constell[i];
    }
    cout<<"constell.txt"<<"矩阵导入完成"<<endl;
    // for (int i = 0; i < row2; i++)
    // {
    //     cout<<"输出H的第"<<i+1<<"行数据"<<endl;
    //     for (int j = 0; j < maxweight1; j++)
    //     {
    //         cout<<n2m[i][j]<<endl;
    //     }
    // }
    //验证一下H*c是不是得0
    flag = gfmatrixmul(H,c,186,372);
    // for (int i = 0; i < row1; i++)
    // {
    //     cout<<flag[i][0]<<endl;
    // }
    
    /*初始化最大似然比函数*/
    float ** Lch = new float *[col1];
    for (int i = 0; i < col1; i++)
    {
        Lch[i] = new float[4];
    }
    int pskdict[8] = {-1,-1,-1,1,1,-1,1,1};
    for (int n = 0; n < col1; n++)
    {
        subconstell = floatslice(constell,n*2,n*2+1);
        for (int i = 0; i<4; i++)
        {
            Lch[n][i] = LLRV(subconstell,pskdict,i);
        }
        delete []subconstell; //一定要释放new出来的内存!!!
    }
    /*初始化Ln2m,Lm2n,Ln2mbuff*/
    float ** Ln2m = new float *[row2];
    float ** Lm2n = new float *[row1];
    float ** Ln2mbuff = new float *[row2];
    for (int i = 0; i < row2; i++)
    {
        Ln2m[i] = new float[maxweight1];
        Ln2mbuff[i] = new float[maxweight1];
    }
    for (int i = 0; i < row1; i++)
    {
        Lm2n[i] = new float[maxweight2];
    }
    return 0;
}

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
void readdata(int row,int col, int **matrix, char *filename)
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
    cout<<filename<<"矩阵导入完成"<<endl;
}
float LLRV(float *subconstell, int *pskdict, int gf)
/*对输入的星座点计算对应的gf元素的对数似然比函数*/
{
    int llrv = 0;
    int index_start = gf*2;
    int index_end = gf*2+1;
    int x_index = 0;
    for (int i = index_start; i <= index_end; i++)
    {
        if(pskdict[i]==1)
        {
            llrv = llrv + subconstell[x_index];
        }
        x_index = x_index + 1;
    }
    return llrv;
}
float *floatslice(float *arry, int start_index, int end_index)
/*用于float类型的切片*/
{
    float *a = new float[end_index-start_index];
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