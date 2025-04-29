#include <iostream>
#include "def.h"
#include "funcs.h"

using namespace std;

/// LB_q  Compute the Lower Bound (LB) for a query
/// This function calculates the lower bound for Dynamic Time Warping (DTW) in nlog(m).
///
/// Variable Explanation:
/// q               : The query sequence.
/// t               : The time series data to compare with the query.
/// m_k             : the length of t.
/// m               : The length of q.
/// lb              : the lower bound (LB).
/// miu_            : Mean values of t.
/// si_             : Standard deviations of t.
/// XX              : Temporary storage for FFT.
/// YY              : Temporary storage for FFT.
/// ZZ              : Temporary storage for FFT.
/// p_forward       : forward plan for performing the FFT transformation.
/// p_backward      : backward plan for performing the inverse FFT transformation.
/// YY_M            : An array for storing the data after M has been padded and then Fourier transformed.
/// YY_UVM          : An array for storing the data after UM_plus_VM has been padded and then Fourier transformed.
/// TM_temp         : Temporary storage for FFT.
/// TU_plus_VM_temp : Temporary storage for FFT.
/// TTM_temp        : Temporary storage for FFT.
/// M_sum           : The sum of all elements in M.
/// c12, c34        : constants, we have already calculated them in main.cpp.
/// special_shared_vector :  A shared vector used to store the calculated distances for the left and right sides(lbkim).
/// TT              : An array for storing the element-wise multiplication of t with itself.
/// D_uv            : ED of UM and VM.
int LB_q(const double q[], double t[], int m_k, int m, double lb[], double threshold, double miu_[], double si_[],
         fftw_complex (*XX), fftw_complex (*YY), fftw_complex (*ZZ), fftw_plan p_forward, fftw_plan p_backward,
         fftw_complex (*YY_M), fftw_complex (*YY_UVM), double *TM_temp, double *TU_plus_VM_temp, double *TTM_temp,
         double M_sum, double c12, double c34, double **special_shared_vector, double TT[], double D_uv)
{
    int i = 0;
    int k_m_l = m_k - m + 1;

    ele_multi(t,t,TT,m_k);

    /// Perform FFT on t and store the result in XX.
    double *TM;
    FFT_x(t, XX, p_forward);

    /// Perform IFFT on the convolution between XX and YY_M and store the result in TM_temp.
    cov_IFFT(XX, YY_M, ZZ, m_k, TM_temp, p_backward);
    /// Adjust TM for the m-th position.
    TM = TM_temp + (m - 1);

    cov_IFFT(XX, YY_UVM, ZZ, m_k, TU_plus_VM_temp, p_backward);
    double *TU_plus_VM = TU_plus_VM_temp + (m - 1);

    FFT_x(TT, YY, p_forward);

    cov_IFFT(YY, YY_M, ZZ, m_k, TTM_temp, p_backward);
    double * TTM = TTM_temp + (m - 1);

    double D_ux2_plus_vx2;
    double threshold2 = threshold * threshold;
    double del = ((2 * threshold + D_uv) * (2 * threshold + D_uv) + D_uv * D_uv) / 2;

    int num_thres = m_k - m + 1;
    double temp1, temp2;

    for (i = 0; i < k_m_l; i++)
    {
        /// Skip subsequences that have an invalid or zero standard deviation.
        if (si_[i] == 0 || isnan(si_[i]))
        {
            lb[i] = threshold + 100;
            continue;
        }
        temp1 = (1 / (si_[i] * si_[i])) * (TTM[i] - 2 * miu_[i] * TM[i] + miu_[i] * miu_[i] * M_sum);
        temp2 = 2 / si_[i];

        D_ux2_plus_vx2 = c12 - temp2 * (TU_plus_VM[i] - miu_[i] * c34) + 2 * temp1;

        /// Perform early abandonment if the distance exceeds the precomputed delta.
        if (D_ux2_plus_vx2 > del) {
            lb[i] = threshold + 100;
            continue;
        }

        /// If the distance is very small, set the lower bound to 0.
        if (D_ux2_plus_vx2 <= D_uv * D_uv + 0.00001)
            lb[i] = 0;
        else {
            double aaa = 2 * D_ux2_plus_vx2 - D_uv * D_uv;
            lb[i] = 0.5 * (sqrt(aaa) - D_uv);
        }

        double d;
        double lb2 = lb[i] * lb[i];

        double dleft;
        double dright;

        double x0 = (t[i] - miu_[i]) / si_[i];
        double y0 = (t[(m - 1 + i)] - miu_[i]) / si_[i];
        dleft = DIST(x0, q[0]);
        dright = DIST(y0, q[m - 1]);

        if (lb2 + dleft + dright >= threshold2) {
            lb[i] = threshold + 100;
            continue;
        }

        /// Continue computing distances for the next points in the subsequence.
        double x1 = (t[(i + 1)] - miu_[i]) / si_[i];
        d = MIN(DIST(x1, q[0]), DIST(x0, q[1]));
        d = MIN(d, DIST(x1, q[1]));
        dleft += d;

        double y1 = (t[(m - 2 + i)] - miu_[i]) / si_[i];
        d = MIN(DIST(y1, q[m - 1]), DIST(y0, q[m - 2]));
        d = MIN(d, DIST(y1, q[m - 2]));
        dright += d;

        if (lb2 + dleft + dright >= threshold2) {
            lb[i] = threshold + 100;
            continue;
        }

        double x2 = (t[(i + 2)] - miu_[i]) / si_[i];
        d = MIN(DIST(x0, q[2]), DIST(x1, q[2]));
        d = MIN(d, DIST(x2, q[2]));
        d = MIN(d, DIST(x2, q[1]));
        d = MIN(d, DIST(x2, q[0]));
        dleft += d;

        double y2 = (t[(m - 3 + i)] - miu_[i]) / si_[i];
        d = MIN(DIST(y0, q[m - 3]), DIST(y1, q[m - 3]));
        d = MIN(d, DIST(y2, q[m - 3]));
        d = MIN(d, DIST(y2, q[m - 2]));
        d = MIN(d, DIST(y2, q[m - 1]));
        dright += d;


        if (lb2 + dleft + dright >= threshold2) {
            lb[i] = threshold + 100;
            continue;
        }

        special_shared_vector[i][0] = dleft;
        special_shared_vector[i][1] = dright;

        /// Decrease the number of subsequences that have passed the threshold check.
        num_thres--;
    }

    /// Return the number of subsequences that pass the threshold check.
    return num_thres;
}


int LB_q_fast(double q[], double t[], int m_k, int m, double lb[], double threshold, double miu_[], double si_[],
              fftw_complex (*Z), fftw_complex (*XX), fftw_complex (*ZZ), fftw_plan p_backward, double UM[], double VM[],
              fftw_complex (*YY_UVM), double *TM, double *TU_plus_VM_temp, double *TTM, double M_module_2, double c12,
              double c34, double **special_shared_vector, double TT[])

{
    int i = 0;
    double  thr_ret=threshold+100;
    int k_m_l = m_k - m + 1;

    double sum_TM = 0;
    double sum_TTM = 0;
    for(i = 0;i < m - 6;i++)
    {
        sum_TM += (t+3)[i];
        sum_TTM += (TT+3)[i];
    }

    double ll = m_k - m +1;
    for(i = 0;i < ll;i++)
    {
        TM[i] = sum_TM;
        TTM[i] = sum_TTM;

        sum_TM = sum_TM - (t+3)[i] + (t+3)[i+m-6];
        sum_TTM = sum_TTM - (TT+3)[i] + (TT+3)[i+m-6];
    }

    cov_IFFT(XX, YY_UVM, ZZ, m_k, TU_plus_VM_temp, p_backward);
    double *TU_plus_VM = TU_plus_VM_temp + (m - 1);

    double threshold2=threshold*threshold;
    double D_uv = ED(UM,VM,m); 
    double del= ((2*threshold+D_uv)*(2*threshold+D_uv) + D_uv*D_uv)/2;
    double D_ux2_plus_vx2;

    double result = 0;
    int num_thres = 0;
    double temp1,temp2,temp3,aaa;
    for(i = 0;i < k_m_l;i++)
    {

        if(lb[i]>threshold) { 
            continue;
        }

        temp1 = ( 1/(si_[i]*si_[i]))*(TTM[i] - 2*miu_[i]*TM[i] + miu_[i]*miu_[i]*M_module_2);
        temp2 = 2/si_[i];
        D_ux2_plus_vx2 = c12- temp2 * (TU_plus_VM[i] - miu_[i]*(c34)) + 2*temp1;

        if (D_ux2_plus_vx2 > del){ 
            lb[i]=thr_ret;
            num_thres ++;
            continue;
        }

        if(D_ux2_plus_vx2<=D_uv*D_uv+0.00001)
            lb[i] = 0;
        else{
            aaa = 2*D_ux2_plus_vx2- D_uv*D_uv;
            lb[i]=0.5 *(sqrt(aaa)- D_uv);

        }

        if(lb[i]*lb[i]+special_shared_vector[i][0]+special_shared_vector[i][1] > threshold2){
            lb[i]=thr_ret; 
            num_thres ++;
        }
    }

    return num_thres;
}

