

#include "funcs.h"
#include "fstream"

using namespace std;

///calculate the tab_q used in  LB_KK_norm_X_on_the_fly_enhance in advance,including some approximate distances
void table_of_q(const double q[], int len, int r, double *result, double Q_min, double derta)
{
    
    double x[21];
    for(int i = 0;i < 21;i++)
    {
        x[i] = Q_min + i*derta;
    }
    double** q_distance = (double**)malloc(20 * sizeof(double*));
    for (int i = 0; i < 20; i++) {
        q_distance[i] = (double*)malloc(len * sizeof(double));
    }
    for(int i = 0;i < 20;i++)
    {
        for(int j = 0;j < len;j++)
        {
            if(q[j] > x[i+1])
                q_distance[i][j] = (x[i+1] - q[j])*(x[i+1] - q[j]);
            else if(q[j] < x[i])
                q_distance[i][j] = (x[i] - q[j])*(x[i] - q[j]);
            else
                q_distance[i][j] = 0;
        }
    }
    
    for(int i = 0;i < 20;i++)
    {
        
        getmin(q_distance[i],len,r,&(result[i*len]));
    }

    for (int i = 0; i < 20; i++) {
        free(q_distance[i]); 
    }
    free(q_distance); 

}