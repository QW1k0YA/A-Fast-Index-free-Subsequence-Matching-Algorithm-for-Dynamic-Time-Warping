#include "def.h"
#include "funcs.h"
typedef double DOUBLE;
using namespace std;
int comp(const void *a, const void* b)
{   Index* x = (Index*)a;
    Index* y = (Index*)b;
    return abs(y->value) - abs(x->value);   
}
int comp_(const void *a, const void* b)
{   Index* x = (Index*)a;
    Index* y = (Index*)b;
    return y->value - x->value;   
}

/// LB_KK_norm_UL_on_the_fly:
///
/// Variable Explanation,
/// q     : the query
/// U, L  : upper and lower envelops for t.
/// seqlen: the length of q
/// miu  : Mean values of t.
/// si   : Standard deviations of t.
/// t     : a circular array keeping the current data.
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
/// order : sorted indices for the query.
/// lbk2  : the final value of LB_KK_norm_UL_on_the_fly
/// special_shared_vector :  A shared vector used to store the calculated distances for the left and right sides(lbkim).
bool LB_KK_norm_UL_on_the_fly(const double q[], const double U[], const double L[], long long seqlen, double miu, double si
,double threshold_2,const int order[],double cb[],double &lb_k2,double special_shared_vector[]) {

    DOUBLE dist = 0.0;
    double temp_u;
    double temp_l;
    
    dist=special_shared_vector[0]+special_shared_vector[1];
    cb[0]=special_shared_vector[0];
    cb[1]=0;
    cb[2]=0;
    cb[seqlen-3]=special_shared_vector[1];
    cb[seqlen-2]=0;
    cb[seqlen-1]=0;

    for (size_t i = 6; i < seqlen; i++) {
        double d = 0;
        temp_u = (U[order[i]] - miu)/si;
        temp_l = (L[order[i]] - miu)/si;
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        dist += d;
        if(threshold_2 < dist)
        {
        
            return 0;

        }
        cb[order[i]] = d;
    }

    lb_k2 = dist;
   
    return 1;
}


/// LB_KK_norm_X_on_the_fly_enhance:
///
/// Variable Explanation,
/// x     : t
/// U, L  : upper and lower envelops for q, which already sorted.
/// seqlen: the length of q
/// miu_  : Mean values of t.
/// si_   : Standard deviations of t.
/// t     : a circular array keeping the current data.
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
/// order : sorted indices for the query.
/// lbk  : the final value of LB_KK_norm_UL_on_the_fly
/// tab_q : the DIST has been already calculated,stored in tab_q
/// special_shared_vector :  A shared vector used to store the calculated distances for the left and right sides(lbkim).
/// Q_min : the MIN value of Q
/// inv_derta : a constant which has been calculated in main.cpp
bool LB_KK_norm_X_on_the_fly_enhance(const double x[], const double U[], const double L[], long long seqlen, double miu,
                                     double si, double threshold_2, const long long int order[], double cb[],
                                     double &lb_k, double *tab_q, double special_shared_vector[], double Q_min,
                                     double inv_derta) {

    DOUBLE dist = 0.0;
    double temp;

    dist=special_shared_vector[0]+special_shared_vector[1];
    cb[0]=special_shared_vector[0];
    cb[1]=0;
    cb[2]=0;
    cb[seqlen-3]=special_shared_vector[1];
    cb[seqlen-2]=0;
    cb[seqlen-1]=0;

    long long ppp;
    for (size_t i = 6; i < seqlen; i++){
        double d = 0;

        temp = (x[order[i]] - miu)/si;

        if (temp >= U[i]) {
            d = (temp - U[i]) * (temp - U[i]);
        } else if (temp <= L[i]) {
            d = (temp - L[i]) * (temp - L[i]);
        }else
        { 
            ppp= int((temp-Q_min)* inv_derta);

            d= tab_q[seqlen* ppp + order[i]];

        }

        dist += d;
        if(threshold_2 < dist)
        {
            return 0;

        }
        cb[order[i]] = d;

    }

    lb_k = dist;
    return 1;
}


DOUBLE  LB_KIM(double t[], double q[],double miu_[],double si_[],double threshold, int m, double special_shared_vector[]){
    int i;
    double  d, dleft,dright;
    double threshold2=threshold*threshold;
    double x0 = (t[i] - miu_[i]) / si_[i];
    double y0 = (t[(m-1+i)] - miu_[i]) / si_[i];
    dleft=DIST(x0, q[0]);
    dright=DIST(y0, q[m - 1]);

    double x1 = (t[(i+1)] - miu_[i]) / si_[i];
    d = MIN(DIST(x1, q[0]), DIST(x0, q[1]));
    d = MIN(d, DIST(x1, q[1]));
    dleft+=d;

    double y1 = (t[(m-2+i)] -  miu_[i]) / si_[i];
    d = MIN(DIST(y1, q[m - 1]), DIST(y0, q[m - 2]) );
    d = MIN(d, DIST(y1, q[m - 2]));
    dright+=d;

    if (dleft+dright  >=threshold2){
        return threshold+100;
    }
    else{

        double x2 = (t[(i+2)] -  miu_[i]) / si_[i];
        d = MIN(DIST(x0, q[2]), DIST(x1, q[2]));
        d = MIN(d, DIST(x2, q[2]));
        d = MIN(d, DIST(x2, q[1]));
        d = MIN(d, DIST(x2, q[0]));
        dleft += d;

        double y2 = (t[(m-3+i)] -  miu_[i]) / si_[i];
        d = MIN(DIST(y0, q[m - 3]), DIST(y1, q[m - 3]));
        d = MIN(d, DIST(y2, q[m - 3]));
        d = MIN(d, DIST(y2, q[m - 2]));
        d = MIN(d, DIST(y2, q[m - 1]));
        dright += d;

        if (dleft+dright  >=threshold2){
            return dleft + dright;
        }
        special_shared_vector[0]=dleft; special_shared_vector[1]=dright;
    }

    return dleft+dright;
}
