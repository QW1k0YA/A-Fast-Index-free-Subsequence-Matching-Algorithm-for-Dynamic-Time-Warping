
#include "def.h"
#include "funcs.h"
#include "string.h"
#include "fstream"

using namespace std;
/// similar to lbq
/// we should calculate the Mask vector in MASK_T
int
LB_t_new(int m_k, int m, double lb[], double u[], double l[], double threshold, double miu_[], double si_[],
         fftw_complex (*XX), fftw_complex (*Y), fftw_complex (*YY), fftw_complex (*ZZ), fftw_plan p_forward,
         fftw_plan p_backward, double *QM_temp, double *UVQM_temp, double *QQM_temp, double *re1, double *re2,
         double *re3, double *M_sum, double *M, double *UM, double *VM, double *UM_VM_plus, double *UM_VM_UM_VM,
         double *UMUM_VMVM, double Q_min, double **special_shared_vector, double lbq[], int bin[],
         double bin_num_divide_T_max_T_min)

{
    m=m-6;
    int k_m_l = m_k - m + 1;

    double  thr_ret=threshold+100;
    mask_T(M, m, m_k, miu_, si_, u, l, Q_min, bin, bin_num_divide_T_max_T_min);

    for(int i = 0;i < m_k;i++)
    {
        UM[i] = u[i]*M[i];
        VM[i] = l[i]*M[i];
        UM_VM_plus[i]=UM[i]+VM[i];
        UM_VM_UM_VM[i] = (UM[i]-VM[i])*(UM[i]-VM[i]);
        UMUM_VMVM[i] = UM[i]*UM[i]+VM[i]*VM[i];
    }

    double *QM;
    double *QQM;
    double *UVQM;

    FFT_x(M, XX, p_forward);
    cov_IFFT(XX, Y, ZZ, m_k, QM_temp, p_backward);
    QM = QM_temp + m - 1;

    cov_IFFT(XX, YY, ZZ, m_k, QQM_temp, p_backward);
    QQM = QQM_temp + m - 1;

    FFT_x(UM_VM_plus, XX, p_forward);
    cov_IFFT(XX, Y, ZZ, m_k, UVQM_temp, p_backward);
    UVQM = UVQM_temp + m - 1;

    mvsum(UM_VM_UM_VM,m_k,m,re1);
    mvsum(UMUM_VMVM,m_k,m,re2);
    mvsum(UM_VM_plus,m_k,m,re3);
    mvsum(M,m_k,m,M_sum);


    double threshold2=threshold*threshold;

    int num_thre = 0;
    for(int i = 3;i < k_m_l-3;i ++)
    {

        int pos=i-3; 
        
        if  (lbq[pos]>threshold) { 
            continue;
        }

        double D_uv = sqrt(re1[i])/si_[pos];
        double D_uvx2 = (1/(si_[pos]*si_[pos]))*(re2[i] - 2*miu_[pos]*re3[i] + 2*miu_[pos]*miu_[pos]*M_sum[i])
                        - (2/si_[pos])*(UVQM[i] - 2 *miu_[pos]*QM[i]) + 2*QQM[i];

        if(D_uvx2<=D_uv*D_uv+0.00001)
            lb[pos]= 0;
        else{
            double aaa=2*D_uvx2- D_uv*D_uv;
            lb[pos]=0.5 *(sqrt(aaa)- D_uv);
        }

        if(lb[pos]*lb[pos]+special_shared_vector[pos][0]+special_shared_vector[pos][1] > threshold2){
            lbq[pos]=thr_ret; 
            num_thre ++;
        }
        
    }

    return num_thre;
}

