#include "llrv.h"
#include "gfcalu.h"
#include "auxil.h"
#include <algorithm>

using namespace std;

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

void Lpostupdate(float *Lpostv, float *Lm2nv, int target)
{
    for (int k = 0; k < 4; k ++)
    {
        Lpostv[k] = Lpostv[k] + Lm2nv[target*4 + k];
    }
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

float *boxplus(float *Ltheta, float *Lro, int H, int q)
{
    int gf0[q-1] = {0};
    int set1_len;
    int *set1;
    float Lpart1 = 0.0;
    float Lpart2 = 0.0;
    float *LL = new float[q];
    for (int k = 0; k < q-1; k ++)
    {
        gf0[k] = k + 1;
        float L4 = Ltheta[k + 1] + Lro[k + 1];
        Lpart2 = max(Lpart2, L4);
    }
    for (int i = 0; i <  q; i ++)
    {
        float L1 = Ltheta[gfdiv(i,H)];
        float L2 = Lro[gfdiv(i,H)];
        Lpart1 = max(L1,L2);
        if (gfdiv(i,H) == 0)
        {
            set1 = gf0;
            set1_len = q - 1;
        }
        else
        {
            set1 = nodediff(gf0, gfdiv(i,H), q-1);
            set1_len = q - 2;
        }
        for (int k = 0; k < set1_len; k ++)
        {
            int x = set1[k];
            float L3 = Ltheta[x] + Lro[gfadd(gfdiv(i,H),x)];
            Lpart1 = max(Lpart1,L3);
        }
        LL[i] = Lpart1 - Lpart2;
    }
    delete []set1;
    return LL;
}