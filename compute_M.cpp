
#include "def.h"
#include "funcs.h"

using namespace std;

///mask data by the order of q
void MASK_DATA(const int *order_of_q, long long m, double* M_data, long long len)
{
     long long m_len = MIN(2 * len, m);
     for(long long i = 0;i < m_len ;i ++)
     {
         M_data[order_of_q[i]] = 1;
     }

    for(long long i = len;i < m;i ++)
    {
        M_data[order_of_q[i]] = 0;
    }

}