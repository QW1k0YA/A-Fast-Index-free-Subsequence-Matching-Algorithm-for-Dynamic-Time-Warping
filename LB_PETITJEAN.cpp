
#include "def.h"
#include "funcs.h"
#define DIST(x,y) ((x-y)*(x-y))
/// a tighter lower bound
bool lbPetitjean_new(double p[], const double q[], double up[], double lp[], const double uq[], const double lq[],
                     const double x[], const double ux[], const double lx[], int w, int m, double threshold_2,
                     double &lbk2, double *cb, double miu, double si) {

    double d = lbk2;
    int l_index = 0;
    int u_index = 0;
    
    int k = 2 * w + 1;
    int q_u[m];
    int q_l[m];
    int head_l = 0, tail_l = -1;
    int head_u = 0, tail_u = -1;

    memcpy(p, q, sizeof(double) * m);
    int flag = 1;

    for(long long i = 3; i < m - 3 + w; i++)
    {

        if(i < m - 3)
        {
            double ux_i = (ux[i] - miu)/si;
            double lx_i = (lx[i] - miu)/si;
            if(q[i] > ux_i)
                p[i] = ux_i;
            else if(q[i] < lx_i)
                p[i] = lx_i;
        }

        if((i >= w))
        {
            if(flag)
            {
                for (int c = 0; c < k - 1 - w; c++) {
                    while (head_l <= tail_l && p[q_l[tail_l]] >= p[c]) tail_l--;
                    q_l[++tail_l] = c;
                    while (head_u <= tail_u && p[q_u[tail_u]] <= p[c]) tail_u--;
                    q_u[++tail_u] = c;
                }

                if(w < 3)
                {
                    for (int b = w ; b < 3 ; b++) {
                        while (head_l <= tail_l && p[q_l[tail_l]] >= p[b]) tail_l--;
                        q_l[++tail_l] = b;
                        while (q_l[head_l] <= b - k) head_l++;
                        lp[l_index++] = p[q_l[head_l]];

                        while (head_u <= tail_u && p[q_u[tail_u]] <= p[b]) tail_u--;
                        q_u[++tail_u] = b;
                        while (q_u[head_u] <= b - k) head_u++;
                        up[u_index++] = p[q_u[head_u]];
                    }
                }
                flag = 0;
            }

            int ii = i - w;
            
            if(i < m)
            {
                while (head_l <= tail_l && p[q_l[tail_l]] >= p[i]) tail_l--;
                q_l[++tail_l] = i;
                while (q_l[head_l] <= i - k) head_l++;
                lp[l_index++] = p[q_l[head_l]];

                while (head_u <= tail_u && p[q_u[tail_u]] <= p[i]) tail_u--;
                q_u[++tail_u] = i;
                while (q_u[head_u] <= i - k) head_u++;
                up[u_index++] = p[q_u[head_u]];
            }
            else
            {
                while (q_u[head_u] <= i - k) head_u++;
                up[u_index++] = p[q_u[head_u]];
                while (q_l[head_l] <= i - k ) head_l++;
                lp[l_index++] = p[q_l[head_l]];
            }
            if(ii >= 3 && ii < m - 3)
            {
                double local_dist=0;

                if(x[ii] > up[ii] && up[ii] > uq[ii])
                {
                     local_dist= (DIST(x[ii], uq[ii]) - DIST(up[ii], uq[ii]));

                }
                else if(q[ii] < lp[ii] && lp[ii] < lq[ii])
                {
                    local_dist = (DIST(x[ii], lq[ii]) - DIST(lp[ii], lq[ii]));
                }
                else if(x[ii] > up[ii])
                {
                    local_dist = DIST(x[ii], up[ii]);
                }
                else if(x[ii] < lp[ii])
                {
                    local_dist = DIST(x[ii], lp[ii]);
                }
                cb[ii]+=local_dist;
                d+=local_dist;
                if(d > threshold_2)
                {
                    return false;
                }

            }
        }
    }
    lbk2 = d;
    return true;

}

