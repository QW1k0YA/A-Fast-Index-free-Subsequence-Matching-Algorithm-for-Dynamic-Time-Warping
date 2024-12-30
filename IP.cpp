///Compute the element-wise product of A and B
double Ip(const double*v1, const double *v2, int len)
{
    double result = 0;
    for(int i = 0;i < len;i++)
    {
        result += v1[i]*v2[i];
    }

    return result;
}