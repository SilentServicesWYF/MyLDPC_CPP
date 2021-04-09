#ifndef LLRV_H
#define LLRV_H

float LLRV(float *subconstell, int *pskdict, int gf);
void Lpostupdate(float *Lpostv, float *Lm2nv, int target);
void LLRVupdate(float *LLRV1, float *LLRV2, int index, int target);
float *LLRVresort(float *LLRV, int H, int q);
float *iterboxplus(float *L, float *Lcnmk, int H, int q);
float *boxplus(float *Ltheta, float *Lro, int H, int q);

#endif