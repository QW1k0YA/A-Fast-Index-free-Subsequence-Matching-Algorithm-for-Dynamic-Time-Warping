#include "funcs.h"

using namespace std;
///Calculating the upper and lower envelopes using monotonic queues
void lower_upper_lemire(double a[], int n, int r, double l[], double u[])
{
    int l_index = 0;
    int u_index = 0;
    int k = 2*r + 1;
    int q_u[n];
    int q_l[n];
    int head_l = 0, tail_l = -1;
    int head_u = 0, tail_u = -1;
    for (int i = 0; i < k - 1 - r; i++) {
        while (head_l <= tail_l && a[q_l[tail_l]] >= a[i]) tail_l--;
        q_l[++tail_l] = i;
        while (head_u <= tail_u && a[q_u[tail_u]] <= a[i]) tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r ; i < n ; i++) {
        while (head_l <= tail_l && a[q_l[tail_l]] >= a[i]) tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k) head_l++;
        l[l_index++] = a[q_l[head_l]];

        while (head_u <= tail_u && a[q_u[tail_u]] <= a[i]) tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k) head_u++;
        u[u_index++] = a[q_u[head_u]];
    }

    for(int i = n;i < n + r;i++)
    {
        while (q_u[head_u] <= i - k) head_u++;
        u[u_index++] = a[q_u[head_u]];
        while (q_l[head_l] <= i - k ) head_l++;
        l[l_index++] = a[q_l[head_l]];
    }

}
///Calculating the lower envelopes using monotonic queues
void getmin(double a[], int n, int r, double l[])
{
    int l_index = 0;
    int u_index = 0;
    int k = 2*r + 1;

    int q_l[n];
    int head_l = 0, tail_l = -1;
    int head_u = 0, tail_u = -1;
    for (int i = 0; i < k - 1 - r; i++) {
        while (head_l <= tail_l && a[q_l[tail_l]] >= a[i]) tail_l--;
        q_l[++tail_l] = i;

    }

    for (int i = k - 1 - r ; i < n ; i++) {
        while (head_l <= tail_l && a[q_l[tail_l]] >= a[i]) tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k) head_l++;
        l[l_index++] = a[q_l[head_l]];

    }
    for(int i = n;i < n + r;i++)
    {
        while (q_l[head_l] <= i - k ) head_l++;
        l[l_index++] = a[q_l[head_l]];
    }
}
void getmax(int k,int q[],double a[],int n,double u[]) {  
    int head = 0, tail = -1;
    for (int i = 1; i < k; i++) {
        while (head <= tail && a[q[tail]] <= a[i]) tail--;
        q[++tail] = i;
    }
    int u_index = 0;
    for (int i = k; i <= n; i++) {
        while (head <= tail && a[q[tail]] <= a[i]) tail--;
        q[++tail] = i;
        while (q[head] <= i - k) head++;
        u[u_index++] = a[q[head]];
    }

    int u_last = u[u_index];
    for(int i = n - k + 1;i < n;i++)
    {
        u[u_index++] = u_last;
    }
}
