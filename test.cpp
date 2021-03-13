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
void readdata(int row,int col, int **matrix, char const *filename);
float LLRV(float *subconstell, int *pskdict, int gf);
float *floatslice(float *arry, int start_index, int end_index);
int *intslice(int *arry, int start_index, int end_index);
int *findnode(int *arry, int node);
void readvector(int row, int *vector, char const *filename);

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
    int n2m_num[col1] = {0};
    int m2n_num[row1] = {0};
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
    readvector(col1,n2m_num,"n2m_num.txt");
    readvector(row1,m2n_num,"m2n_num.txt");
    ifstream fin;
    fin.open("constell.txt");
    for (int i = 0; i < 2*row2; i++)
    {
        fin >> constell[i];
    }
    cout<<"constell.txt"<<"矩阵导入完成"<<endl;
    // for (int i = 0; i < 372*2; i++)
    // {
    //     cout<<constell[i]<<endl;
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
        // for (int k = 0;  k < 2; k++)
        // {
        //     cout<<subconstell[k]<<endl;
        // }
        for (int i = 0; i<4; i++)
        {
            Lch[n][i] = LLRV(subconstell,pskdict,i);
        }
        delete []subconstell; //一定要释放new出来的内存!!!
    }
    // for (int i = 0; i < col1; i++)
    // {
    //     cout<<"输出Lch的第"<<i+1<<"行数据"<<endl;
    //     for (int j = 0; j < 4; j++)
    //     {
    //         cout<<Lch[i][j]<<endl;
    //     }
    // }
    /*初始化Ln2m,Lm2n,Ln2mbuff*/
    float ** Ln2m = new float *[row2];
    float ** Lm2n = new float *[row1];
    float ** Ln2mbuff = new float *[row2];
    for (int i = 0; i < row2; i++)
    {
        Ln2m[i] = new float[maxweight1];
        Ln2mbuff[i] = new float[maxweight1*4];
    }
    for (int i = 0; i < row1; i++)
    {
        Lm2n[i] = new float[maxweight2];
    }
    for (int k = 0; k < col1; k++)
    {
        int avm_index = 0;
        while ((n2m[k][avm_index] != 0) && (avm_index < maxweight1))
        {
            int Lch_index = 0;
            for (int buff_index = avm_index*4; buff_index < (avm_index + 1)*4; buff_index++)
            {
                Ln2mbuff[k][buff_index] = Lch[k][Lch_index];
                Lch_index = Lch_index + 1;
            }
            avm_index = avm_index + 1;
        }
    }
    // for (int i = 0; i < 10; i++)
    // {
    //     cout<<"输出Ln2mbuff的第"<<i+1<<"行数据"<<endl;
    //     for (int j = 0; j < maxweight1*4; j++)
    //     {
    //         cout<<Ln2mbuff[i][j]<<endl;
    //     }
    // }
    float ** Lpost = new float *[col1]; //为Lpost申请内存
    int ** est_c = new int *[col1]; //为最大似然估计估计出来的码元申请内存
    for (int k = 0; k < col1; k ++)
    {
        Lpost[k] = new float[4];
        est_c[k] = new int[1];
    }
    /* 迭代开始 */
    int maxiter = 5000;
    for (int iter = 0; iter < maxiter; iter ++)
    {
        // 尝试性解码
        //为Lpost赋值
        for (int k = 0; k < col1; k ++)
        {
            for (int j = 0; j < 4; j ++)
            {
                Lpost[k][j] = Lch[k][j];
            }
        }
    //更新Lpost的值
        for (int n = 0; n < 372; n ++)
        {
            
        }
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
    cout<<filename<<"矩阵导入完成"<<endl;
}
void readvector(int row, int *vector, char const *filename)
{
    ifstream fin;
    fin.open(filename);
    for (int i = 0; i < row; i++)
    {
        fin >> vector[i];
    }
    cout<<filename<<"向量导入完成"<<endl;
}
float LLRV(float *subconstell, int *pskdict, int gf)
/*对输入的星座点计算对应的gf元素的对数似然比函数*/
{
    float llrv = 0;
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
int *findnode(int **arry, int node)
{
    int *temp;
    temp = arry[node];
    int len = sizeof(temp)/sizeof(temp[0]);
    int count = 0;
    for (int k = 0; k < len; k++)
    {
        if (temp[k] != 0)
        {
            count = count + 1;
        }
    }
    int *targetnode = new int[count];
    for (int k = 0; k < count; k ++)
    {
        targetnode[k] = temp[k];
    }
    return targetnode;
}