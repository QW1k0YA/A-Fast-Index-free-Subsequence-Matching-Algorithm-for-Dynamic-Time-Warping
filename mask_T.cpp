
#include "def.h"
#include "funcs.h"
#include "string.h"
#include "fstream"
using namespace std;

///calculate the index in bin
///The array 'bin' has already been initialized in main.cpp
int H(double v, double bin_num_divide_T_max_T_min, double Q_min, double bin_num)
{
    
    int result = 0;
    result = (v - Q_min) * bin_num_divide_T_max_T_min;

    if(result < 0)
        result = 0;
    if(result > (bin_num-1))
        result = bin_num - 1;
    return result;
}

///calculate the MASK vector used in lbt
void mask_T(double *M, int m, int m_k, double *miu, double *si, double *u, double *v, double Q_min, int bin[],
            double bin_num_divide_T_max_T_min)
{
    
    double bin_num = BIN ;

    double u_p;
    double v_p;

    int L = m;

    int L_2 = L/2;

    for(int i = 0;i < m_k;i++)
    {
        int tt = i - L_2;
        if((i >= (L_2 )) && (i <= (m_k - L + L_2)))
        {
            u_p = (u[i] - miu[tt])/si[tt];
            v_p = (v[i] - miu[tt])/si[tt];
        }
        else if(i < L_2)
        {
            u_p = (u[i] - miu[0])/si[0];
            v_p = (v[i] - miu[0])/si[0];
        }
        else
        {
            u_p = (u[i] - miu[m_k - m])/si[m_k - m];
            v_p = (v[i] - miu[m_k - m])/si[m_k - m];
        }

        if((bin[H(u_p+0.1,bin_num_divide_T_max_T_min,Q_min,bin_num)] - bin[H(v_p-0.1,bin_num_divide_T_max_T_min,Q_min,bin_num)]) < (0.5*L))
        {
            M[i] = 1;
           
        }
        else
        {
            M[i] = 0;
            
        }
    }

}