
#include "funcs.h"

using namespace std;
///compute the dot product of A and B
double dot_mul(const double *A,const double *B,long long len)
{
    double sum = 0;
    for(size_t i = 0;i < len;i++)
    {
        sum += (A[i] * B[i]);
    }
    return sum;
}