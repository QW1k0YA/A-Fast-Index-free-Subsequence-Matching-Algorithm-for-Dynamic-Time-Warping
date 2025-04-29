
#include "def.h"
#include "funcs.h"
#define DIST(x,y) ((x-y)*(x-y))
/// similar to lbq
long long int
lb_data(double *lb, double *M, double *q, long long m, double *U, double *L, double *miu, double *MQ, long long m_k,
        fftw_complex (*V1), fftw_plan p_forward, fftw_plan p_backward, double *u_plus_l, double *U_2, double *L_2,
        double *temp_mk4, fftw_complex (*V2), fftw_complex (*V3), double *b_temp, double *d, double *f, double *si,
        double *miu_2, double *U_mul_L, double *d_temp, double *e_temp, double *f_temp, double threshold2,
        double **special_shared_vector, double *lbq)
{
    
    M[0]=0;M[1]=0;M[2]=0;M[m-1]=0;M[m-2]=0;M[m-3]=0;

    for(long long i = 0;i < m_k - m + 1;i++)
    {
        miu_2[i] = miu[i]*miu[i];
    }
    ele_multi(M,q,MQ,m);
    double a = dot_mul(MQ, MQ, m);

    ele_plus(U, L, u_plus_l, m_k);
    inv_and_padding_y(m_k, m, MQ, V1, p_forward);
    FFT_x(u_plus_l, V2, p_forward);
    cov_IFFT(V2, V1, V3, m_k, b_temp, p_backward);
    double *b = b_temp + m - 1;

    double c_temp = dot_mul(M,q,m);

    ele_multi(U, U, U_2, m_k);
    ele_multi(L, L, L_2, m_k);
    ele_plus(U_2, L_2, temp_mk4, m_k);
    inv_and_padding_y(m_k, m, M, V1, p_forward);

    cov_IFFT(V2, V1, V3, m_k, e_temp, p_backward);
    double *e_ = e_temp + m -1;

    FFT_x(temp_mk4, V2, p_forward);
    cov_IFFT(V2, V1, V3, m_k, d_temp, p_backward);
    d = d_temp + m - 1;

    ele_multi(U, L, U_mul_L, m_k);
    FFT_x(U_mul_L, V2, p_forward);
    cov_IFFT(V2, V1, V3, m_k, f_temp, p_backward);
    f = f_temp + m -1;

    double M_module_2 = 0;
    for(long long i = 0;i < m;i++)
    {
        M_module_2 += M[i];
    }

    double lb1;
    double lb2;

    long long num_thre = 0;

    for(long long i = 0;i < m_k - m + 1;i++)
    {
        if (lbq[i]*lbq[i]>threshold2) {
            continue;
        }
        lb1 = 2*a - (2/si[i])*(b[i] - 2*(miu[i]*c_temp)) + (1/(si[i]*si[i]))*(d[i] - 2*(miu[i]*e_[i])+2*(M_module_2 * miu_2[i]));
        lb2 = (1/(si[i]*si[i])) * (d[i] - 2*f[i]);
        lb[i] = 0.5 * (sqrt(2*lb1 - lb2) - sqrt(lb2));
        
        if(lb[i]*lb[i]+special_shared_vector[i][0]+special_shared_vector[i][1] > threshold2)
        {
            num_thre ++;
        }
        lbq[i] = MAX(lbq[i], lb[i]);
    }

    return num_thre;
}