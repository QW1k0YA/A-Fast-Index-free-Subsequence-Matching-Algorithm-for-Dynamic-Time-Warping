
#ifndef SCALE_DTW_FUNCS_H
#define SCALE_DTW_FUNCS_H
#include "def.h"
#include "vector"
#include "fftw3.h"

extern size_t lb_kill_early_abandon_counter;

void ele_multi(const double *v1,const double *v2,double *re,int len);
double ED(const double* u,const double* v,int len);
double Ip(const double*v1, const double *v2, int len);
void ele_plus(const double v1[],const double v2[],
              double re[],int len);

int H(double v, double bin_num_divide_T_max_T_min, double Q_min, double bin_num);
void mask_T(double *M, int m, int m_k, double *miu, double *si, double *u, double *v, double Q_min, int bin[],
            double bin_num_divide_T_max_T_min);

double normal_cdf(double x);
void mvmean(double *a,int len_a,int l,double *miu,double *si);
void mvsum(double *a,int len_a,int l,double *result);
void  inv_and_padding_y(int n,int m,double *y, fftw_complex * YY,fftw_plan p);
void FFT_x(double *x, fftw_complex (*XX), fftw_plan p);
void cov_IFFT(fftw_complex (*XX), fftw_complex (*YY), fftw_complex (*ZZ), int n, double *z, fftw_plan p);

bool LB_KK_norm_X_on_the_fly_enhance(const double x[], const double U[], const double L[], long long seqlen, double miu,
                                     double si, double threshold_2, const long long int order[], double cb[],
                                     double &lb_k, double *tab_q, double special_shared_vector[], double Q_min,
                                     double inv_derta);
bool LB_KK_norm_UL_on_the_fly(const double q[], const double U[], const double L[], long long seqlen, double miu, double si
        ,double threshold_2,const int order[],double cb[],double &lb_k2,double special_shared_vector[]) ;
double  LB_KIM(double t[], double q[],double miu_[],double si_[],double threshold, int m, double special_shared_vector[]);
double dtw(double* A, double* B,  int m, int r, double threshold_2,double *cb);
void lower_upper_lemire(double a[], int n, int r, double l[], double u[]);

int LB_q(const double q[], double t[], int m_k, int m, double lb[], double threshold, double miu_[], double si_[],
         fftw_complex (*XX), fftw_complex (*YY), fftw_complex (*ZZ), fftw_plan p_forward, fftw_plan p_backward,
         fftw_complex (*YY_M), fftw_complex (*YY_UVM), double *TM_temp, double *TU_plus_VM_temp, double *TTM_temp,
         double M_sum, double c12, double c34, double **special_shared_vector, double TT[], double D_uv);
int LB_q_fast(double q[], double t[], int m_k, int m, double lb[], double threshold, double miu_[], double si_[],
              fftw_complex (*Z), fftw_complex (*XX), fftw_complex (*ZZ), fftw_plan p_backward, double UM[], double VM[],
              fftw_complex (*YY_UVM), double *TM, double *TU_plus_VM_temp, double *TTM, double M_module_2, double c12,
              double c34, double **special_shared_vector, double TT[]);
int
LB_t_new(int m_k, int m, double lb[], double u[], double l[], double threshold, double miu_[], double si_[],
         fftw_complex (*XX), fftw_complex (*Y), fftw_complex (*YY), fftw_complex (*ZZ), fftw_plan p_forward,
         fftw_plan p_backward, double *QM_temp, double *UVQM_temp, double *QQM_temp, double *re1, double *re2,
         double *re3, double *M_sum, double *M, double *UM, double *VM, double *UM_VM_plus, double *UM_VM_UM_VM,
         double *UMUM_VMVM, double Q_min, double **special_shared_vector, double lbq[], int bin[],
         double bin_num_divide_T_max_T_min);

void getmin(double a[], int n, int r, double l[]);
void table_of_q(const double q[], int len, int r, double *result, double Q_min, double derta);
int comp(const void *a, const void* b);
int comp_(const void *a, const void* b);
double dot_mul(const double *A,const double *B,long long len);
void MASK_DATA(const int *order_of_q, long long m, double* M_data, long long len);
long long int
lb_data(double *lb, double *M, double *q, long long m, double *U, double *L, double *miu, double *MQ, long long m_k,
        fftw_complex (*V1), fftw_plan p_forward, fftw_plan p_backward, double *u_plus_l, double *U_2, double *L_2,
        double *temp_mk4, fftw_complex (*V2), fftw_complex (*V3), double *b_temp, double *d, double *f, double *si,
        double *miu_2, double *U_mul_L, double *d_temp, double *e_temp, double *f_temp, double threshold2,
        double **special_shared_vector, double *lbq);
bool lbPetitjean_new(double p[], const double q[], double up[], double lp[], const double uq[], const double lq[],
                     const double x[], const double ux[], const double lx[], int w, int m, double threshold_2,
                     double &lbk2, double *cb, double miu, double si);
double MON_dtw(
        const double* lines,
        const double* cols,
        const double* cb,
        int l,
        int w,
        double bsf
);

#endif 
