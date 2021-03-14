#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>

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
int findnode(int **arry, int node1, int node2,int len);
void readvector(int row, int *vector, char const *filename);
void Lpostupdate(float *Lpostv, float *Lm2nv, int target);
int *nodediff(int *mset, int m, int n2m_num);
void LLRVupdate(float *LLRV1, float *LLRV2, int index, int target);
float *LLRVresort(float *LLRV, int H, int q);
void boxplusupdate(float *L, float *Lcmnk, int H, int q);
float *iterboxplus(float *L, float *Lcnmk, int H, int q);

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
    delete []flag;
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
    int iterflag = 1;
    int iter_num = 0;
    int maxiter = 5000;
    while(iterflag == 1 && iter_num < maxiter)
    {
        iter_num = iter_num + 1;
        /*尝试性解码*/
        //为Lpost赋值
        for (int k = 0; k < col1; k ++)
        {
            for (int j = 0; j < 4; j ++)
            {
                Lpost[k][j] = Lch[k][j];
                // cout<<Lpost[k][j]<<endl;
            }
        }
        //更新Lpost的值
        for (int n = 0; n < 372; n ++)
        {
            //搜寻所有与当前n连接的节点的m节点
            int mset[n2m_num[n]] = {0};
            for (int k = 0; k < n2m_num[n]; k ++)
            {
                mset[k] = n2m[n][k];
            }
            //搜索到所有与n连接的m节点之后对每一个被连接的m节点搜索其连接的n节点在Lm2n中的位置然后更新Lpost
            for (int k = 0; k < n2m_num[n]; k ++) //n2m_num[n]就是当前被连接的m节点的数量
            {
                //搜寻被当前m节点连接的n在Lm2n中的位置
                int targetm2n = findnode(m2n,mset[k]-1,(n+1),5); //第三个参数要+1是因为此参数是被搜寻的n而不是角标了
                //而第二个参数则是以角标形式出现的
                Lpostupdate(Lpost[n],Lm2n[mset[k]-1],targetm2n);
                // for (int j = 0;j < 4;j ++)
                // {
                //     cout<<Lpost[n][j]<<endl;
                // }
            }
            est_c[n][0] = max_element(Lpost[n],Lpost[n]+4) - Lpost[n];  //根据最大后验估计技术est_c
        }
        // for (int i = 0; i < col1; i++)
        // {
        //     cout<<"输出Lpost的第"<<i+1<<"行数据"<<endl;
        //     for (int j = 0; j < 4; j++)
        //     {
        //         cout<<Lpost[i][j]<<endl;
        //     }
        // }
        //统计错误码元数
        int errorcount = 0;
        for (int k = 0; k < col1; k ++)
        {
            // cout<<(est_c[k][0]-c[k][1])<<endl;
            if ((est_c[k][0]-c[k][0]) != 0)
            {
                errorcount = errorcount + 1;
            }
        }
        cout<<"错误码元数："<<errorcount<<endl;
        //判断是否解码成功条件
        flag = gfmatrixmul(H,est_c,186,372);
        for (int k = 0; k < 186; k ++)
        {
            if (flag[k] != 0)
            {
                iterflag = 1;
                break;
            }
            iterflag = 0;
        }
        delete []flag;
        //水平信息传递
        //初始化Ln2m
        for (int n = 0; n < col1; n ++)
        {
            for (int k = 0; k < maxweight1*4; k ++)
            {
                Ln2m[n][k] = Ln2mbuff[n][k];
            }
        }
        for (int n = 0; n < col1; n ++)
        {
            //搜寻所有与当前n连接的节点的m节点
            int mset[n2m_num[n]] = {0};
            for (int k = 0; k < n2m_num[n]; k ++)
            {
                mset[k] = n2m[n][k];
            }
            //对每个n连接的m节点计算除了本m以外连接到n的m节点进行LLRV更新
            for (int k = 0; k < n2m_num[n]; k ++)
            {
                int *msubset;
                msubset = nodediff(mset, mset[k], n2m_num[n]);
                for (int mm = 0; mm < n2m_num[n] - 1; mm ++)
                {
                    int targetm2n = findnode(m2n,msubset[mm]-1,(n + 1),5); //元素转角标必须-1，角标转元素必须+!
                    LLRVupdate(Ln2m[n], Lm2n[n], k, targetm2n);
                }
            }  
        }
        // 垂直信息传递
        for (int m = 0; m < row1; m ++)
        {
            //找出与本n相连的m节点
            int nset[m2n_num[m]];
            for (int k = 0; k < m2n_num[m]; k ++)
            {
                nset[k] = m2n[m][k];
            }
            // 初始化Ltheta Lro 后面别忘删
            float ** Ltheta = new float *[m2n_num[m]]; //Ltheta最后一行无意义
            float ** Lro = new float *[m2n_num[m]]; //Lro第一行无意义
            // 初始化Ltheta[0]
            int targetn2m = findnode(n2m, nset[0], (m + 1), 3); //寻找当前m在Ln2m第n行中的角标
            float *Lcnm1 = floatslice(Ln2m[nset[0] - 1], targetn2m*4, (targetn2m + 1)*4);
            Ltheta[0] = LLRVresort(Lcnm1, H[m][nset[0] - 1], 4); // 节点转角标-1,内存在Ltheta结束后释放
            delete []Lcnm1; //Lcnm1在LLRVresort就可以释放了
            // 初始化Lro[0]
            targetn2m = findnode(n2m, nset[m2n_num[m]], (m + 1), 3);
            float *Lcnmend = floatslice(Ln2m[nset[m2n_num[m] - 1]], targetn2m*4, (targetn2m + 1)*4);
            Lro[m2n_num[m] - 1] = LLRVresort(Lcnmend, H[m][nset[m2n_num[m] - 1] - 1], 4);
            delete []Lcnmend;
            // 为Ltheta Lro循环赋值
            for (int k = 1; k < m2n_num[m]; k ++)
            {
                int target1 = findnode(n2m, nset[k], (m + 1), 3); //Ltheta的target
                int target2 = findnode(n2m, nset[m2n_num[m] -1 - k] ,(m + 1), 3); //Lro的target
                float *Lcnmk = floatslice(Ln2m[nset[k] - 1], target1*4, (target1 + 1)*4);
                Ltheta[k] = iterboxplus(Ltheta[k-1], Lcnmk, H[m][nset[k]-1], 4);
                Lcnmk = floatslice(Ln2m[nset[m2n_num[m] -1 - k] - 1], target2*4, (target2 + 1)*4);
                Lro[m2n_num[m] - 1 - k] = iterboxplus(Lro[m2n_num[m] - k], Lcnmk, H[m][nset[m2n_num[m] -1 - k]-1], 4);
                delete []Lcnmk;
            }
            for (int k = 0; k < m2n_num[m]; k ++)
            {
                float *Lm2ntemp;
                if (k == 0)
                {
                    Lm2ntemp = LLRVresort(Lro[k + 1], gfdiv(1,H[m][nset[k] - 1]), 4);
                }
                else if (k == m2n_num[m] - 1)
                {
                    Lm2ntemp = LLRVresort(Ltheta[k - 1],gfdiv(1,H[m][nset[k] - 1]), 4);
                }
                else
                {
                    // Lm2ntemp = 
                }
            }

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
void Lpostupdate(float *Lpostv, float *Lm2nv, int target)
{
    for (int k = 0; k < 4; k++)
    {
        Lpostv[k] = Lpostv[k] + Lm2nv[target*4 + k];
    }
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
void LLRVupdate(float *LLRV1, float *LLRV2, int index, int target)
{
/*用于更新Ln2m*/
    for (int k = 0; k < 4; k ++)
    {
        LLRV1[index*4 + k] = LLRV1[index*4 + k] + LLRV2[target*4 + k];
    }
}
float *LLRVresort(float *LLRV, int H, int q)
{
    float *templlrv = new float[q];
    templlrv[0] = LLRV[0]; //0除H永远得0
    for (int k = 1; k < q; k ++)
    {
        int gf_index = gfdiv(k,H);
        templlrv[k] = LLRV[gf_index]; 
    }
    return templlrv;
}
float *iterboxplus(float *L, float *Lcnmk, int H, int q)
/*迭代计算Ltheta,Lro*/
{
    int gf0[q-1] = {0};
    int set1_len;
    int *set1;
    float Lpart1 = 0.0;
    float Lpart2 = 0.0;
    float *LL = new float[q]; //新变量
    for (int k = 0; k < q-1; k++)
    {
        gf0[k] = k + 1;
        float L4 = L[k + 1] + Lcnmk[gfdiv(k + 1,H)];
        Lpart2 = max(Lpart2,L4);
    }
    for (int i = 0; i < q; i ++)
    {
        float L1 = L[i];
        float L2 = Lcnmk[gfdiv(i,H)];
        Lpart1 = max(L1,L2);
        if (i == 0)
        {
            set1 = gf0;
            set1_len = q - 1;
        }
        else
        {
            set1 = nodediff(gf0, i, q-1); //set1用完删除
            set1_len = q - 2;
        }
        for (int k = 0; k < set1_len; k ++)
        {
            int x = set1[k];
            float L3 = L[x] + Lcnmk[gfdiv(gfadd(i,x),H)];
            Lpart1 = max(Lpart1,L3);
        }
        LL[i] = Lpart1 - Lpart2;
    }
    delete []set1; //清除set1的内存
    return LL;
}
float *boxplus(float *ltheta, float *Lro, float H, int q)
{
    
}