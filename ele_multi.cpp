///Compute the element-wise product of A and B
void ele_multi(const double *v1,const double *v2,double *re,int len)
{
    for(int i = 0;i < len;i++)
    {
        re[i] = v1[i]*v2[i];
    }

}