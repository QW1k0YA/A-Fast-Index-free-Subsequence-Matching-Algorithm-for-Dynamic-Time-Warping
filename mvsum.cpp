///calculates the moving sum of a sliding window of size l over an array a
void mvsum(double *a,int len_a,int l,double *result)
{
    double sum = 0;
    for(int i = 0;i < l;i++)
    {
        sum+= a[i];
    }

    double ll = len_a - l +1;
    for(int i = 0;i < ll;i++)
    {
        result[i] = sum;
        sum = sum - a[i] + a[i+l];
    }

}
