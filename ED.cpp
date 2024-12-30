
#include "cmath"

double ED(const double* u,const double* v,int len)
{
    double result = 0;
    for(int i = 0;i < len;i++)
    {
        result += (u[i]-v[i])*(u[i]-v[i]);
        
    }
    result = sqrt(result);
    return result;
}

