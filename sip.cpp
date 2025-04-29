
#include "fftw3.h"
#include "funcs.h"

///padding and then FFT
void  inv_and_padding_y(int n,int m,double *y, fftw_complex * YY,fftw_plan p)

{
    double y_temp[n];
    for(int i = 0 ; i < n ; i++ )
    {
        if(i < m )
            y_temp[i] = y[m-i-1]; 
        else
            y_temp[i] = 0;
    }
    fftw_execute_dft_r2c(p,y_temp,YY);
   
}
///Perform FFT on x
void FFT_x(double *x, fftw_complex (*XX), fftw_plan p)
{
    fftw_execute_dft_r2c(p,x,XX);
}
/// Perform IFFT on the convolution between XX and YY,and store the result in z
void cov_IFFT(fftw_complex (*XX), fftw_complex (*YY), fftw_complex (*ZZ), int n, double *z, fftw_plan p)
{
    ///For real-valued sequences, the FFT results are always conjugate symmetric, so only half of the points need to be computed.
    int i_ceil = n/2;
    double inv_n = 1.0/n;

    ///If n is odd, the middle value must also be computed.
    for(int i = 0 ; i <= i_ceil; i++)
    {
        ///computing the real part
        ZZ[i][0] = (XX[i][0]*YY[i][0] - XX[i][1]*YY[i][1])*inv_n;
        ///computing the imaginary part
        ZZ[i][1] = (XX[i][1]*YY[i][0] + XX[i][0]*YY[i][1])*inv_n;
    }

    fftw_execute_dft_c2r(p,ZZ,z);
}

