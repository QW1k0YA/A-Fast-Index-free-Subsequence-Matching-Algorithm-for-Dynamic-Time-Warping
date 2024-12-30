
#include "def.h"

///Calculate the cumulative distribution function (CDF) of the normal distribution
double normal_cdf(double x) {

    if(x >= 0)
    {
        return cdf[int(x*999/5)];
    }
    else
    {
        return 1 - cdf[int((-1)*x*999/5)];
    }
}

