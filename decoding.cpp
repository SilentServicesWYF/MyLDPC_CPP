#include "gfcalu.h"
#include "auxil.h"
#include "llrv.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <time.h>
#include <sys/time.h>

using namespace std;

int main()
{
    //参数设定
    int row1 = 1344;
    int col1 = 2688;
    int row2 = 2688;
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
    float ax_L[4] = {0};
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
    readdata(row1,col1,H,"data/H.txt");
    readdata(row2,col2,c,"data/c.txt");
    readdata(row1,maxweight2,m2n,"data/m2n.txt");
    readdata(row2,maxweight1,n2m,"data/n2m.txt");
    readvector(col1,n2m_num,"data/n2m_num.txt");
    readvector(row1,m2n_num,"data/m2n_num.txt");
    ifstream fin;
    fin.open("data/constell.txt");
    for (int i = 0; i < 2*row2; i++)
    {
        fin >> constell[i];
    }
    std::cout<<"constell.txt"<<"矩阵导入完成"<<endl;
    //验证一下H*c是不是得0
    flag = gfmatrixmul(H,c,row1,col1);
    delete []flag;
    float ** Lch = new float *[col1];
    int pskdict[8] = {-1,-1,-1,1,1,-1,1,1};
    for (int n = 0; n < col1; n++)
    {
        Lch[n] = new float[4];
        subconstell = floatslice(constell,n*2,n*2+1);
        for (int i = 0; i<4; i++)
        {
            Lch[n][i] = LLRV(subconstell,pskdict,i);
        }
        delete []subconstell; //一定要释放new出来的内存!!!
    }
    // for (int i = 0; i < col1; i ++)
    // {
    //     cout<<"d_Lch第"<<i+1<<"行数据"<<endl;
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
        Ln2m[i] = new float[maxweight1*4];
        Ln2mbuff[i] = new float[maxweight1*4];
    }
    for (int i = 0; i < row1; i++)
    {
        Lm2n[i] = new float[maxweight2*4];
    }
    // 为Ln2mbuff赋值
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
    int maxiter = 50000;

    clock_t decodebegin = clock();
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
            }
        }
        //更新Lpost的值
        for (int n = 0; n < col1; n ++)
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
            }
            est_c[n][0] = max_element(Lpost[n],Lpost[n]+4) - Lpost[n];  //根据最大后验估计技术est_c
        }
        // 统计错误码元数
        int errorcount = 0;
        for (int k = 0; k < col1; k ++)
        {
            if ((est_c[k][0]-c[k][0]) != 0)
            {
                errorcount = errorcount + 1;
            }
        }
        cout<<"迭代第"<<iter_num<<"次"<<"错误码元数："<<errorcount<<endl;
        // 判断是否解码成功条件
        flag = gfmatrixmul(H,est_c,row1,col1);
        iterflag = 0;
        for (int k = 0; k < row1; k ++)
        {
            if (flag[k][0] != 0)
            {
                iterflag = 1;
                break;
            }
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
                    int targetm2n = findnode(m2n,msubset[mm]-1,(n + 1),5); //元素转角标必须-1，角标转元素必须+1
                    LLRVupdate(Ln2m[n], Lm2n[msubset[mm] - 1], k, targetm2n);
                }
                delete []msubset;
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
            int targetn2m = findnode(n2m, nset[0] - 1, (m + 1), 3); //寻找当前m在Ln2m第n行中的角标
            float *Lcnm1 = floatslice(Ln2m[nset[0] - 1], targetn2m*4, (targetn2m + 1)*4 - 1);
            Ltheta[0] = LLRVresort(Lcnm1, H[m][nset[0] - 1], 4); // 节点转角标-1,内存在Ltheta结束后释放
            delete []Lcnm1; //Lcnm1在LLRVresort就可以释放了
            // 初始化Lro[0]
            targetn2m = findnode(n2m, nset[m2n_num[m] - 1] - 1, (m + 1), 3);
            float *Lcnmend = floatslice(Ln2m[nset[m2n_num[m] - 1] - 1], targetn2m*4, (targetn2m + 1)*4 - 1);
            Lro[m2n_num[m] - 1] = LLRVresort(Lcnmend, H[m][nset[m2n_num[m] - 1] - 1], 4);
            delete []Lcnmend;
            // 为Ltheta循环赋值
            for (int k = 1; k < m2n_num[m]-1; k ++)
            {
                int target = findnode(n2m, nset[k] - 1, (m + 1), 3); //Ltheta的target
                float *Lcnmk = floatslice(Ln2m[nset[k] - 1], target*4, (target + 1)*4 - 1);
                Ltheta[k] = iterboxplus(Ltheta[k-1], Lcnmk, H[m][nset[k]-1], 4);
                delete []Lcnmk;
            }
            Ltheta[m2n_num[m]-1] = ax_L;
            // 为Lro循环赋值
            for (int k = m2n_num[m] - 2; k > 0; k --)
            {
                int target = findnode(n2m, nset[k] - 1 ,(m + 1), 3); //Lro的target
                float *Lcnmk = floatslice(Ln2m[nset[k] - 1], target*4, (target + 1)*4 - 1);
                Lro[k] = iterboxplus(Lro[k + 1], Lcnmk, H[m][nset[k]-1], 4);
            }
            Lro[0] = ax_L;
            for (int k = 0; k < m2n_num[m]; k ++)
            {
                float *Lm2ntemp;
                if (k == 0)
                {
                    Lm2ntemp = LLRVresort(Lro[1], gfdiv(1,H[m][nset[k] - 1]), 4);
                }
                else if (k == m2n_num[m] - 1)
                {
                    Lm2ntemp = LLRVresort(Ltheta[m2n_num[m] - 2], gfdiv(1,H[m][nset[k] - 1]), 4);
                }
                else
                {
                    Lm2ntemp = boxplus(Ltheta[k - 1], Lro[k + 1], gfdiv(1,H[m][nset[k] - 1]), 4); 
                }
                for (int j = 0; j < 4; j ++)
                {
                    Lm2n[m][k*4 + j] = Lm2ntemp[j];
                }
            }
            delete []Ltheta;
            delete []Lro;
        }
    }
    if (iter_num < maxiter)
    {
        cout<<"解码成功"<<endl;
    }
    clock_t decodeend = clock();
    double timeuse = ((double)(decodeend - decodebegin))/CLOCKS_PER_SEC;
    std::cout<<iter_num<<"轮解码耗时"<<timeuse<<"s"<<endl;
    return 0;
}

