
#include "cmath"
using namespace std;
/// calculates the moving average (mean) and standard deviation of a sliding window over an array
void mvmean(double *a,int len_a,int l,double *miu,double *si)
{

    double sum1 = 0;
    double sum2 = 0;
    for(int i = 0;i < l;i++)
    {
        sum1+= a[i];
        sum2+= a[i]*a[i];
    }

    double ll = len_a - l +1;
    for(int i = 0;i < ll;i++)
    {
        miu[i] = sum1/l;
        si[i] = sqrt(sum2/l-miu[i]*miu[i]);
        sum1 = sum1 - a[i] + a[i+l];
        sum2 = sum2 - a[i]*a[i] + a[i+l]*a[i+l];
    }

}

