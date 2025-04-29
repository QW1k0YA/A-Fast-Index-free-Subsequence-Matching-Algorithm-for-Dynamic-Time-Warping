
#include "def.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
#define POSITIVE_INFINITY 1E20

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double* A, double* B,  int m, int r, double threshold_2,double *cb)
{

    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;

    cost = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=INF;

    cost_prev = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=INF;

    for (i=0; i<m; i++)
    {
        k = MAX(0, r - i);
        min_cost = INF;

        for(j=MAX(0, i - r); j <= MIN(m - 1, i + r); j++, k++)
        {
            
            if ((i==0)&&(j==0))
            {
                cost[k]=DIST(A[0], B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];

            cost[k] = MIN(MIN(x, y) , z) + DIST(A[i], B[j]);

            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            }
        }

        if (i+r < m-1 && min_cost + cb[i+r+1] >= threshold_2)
        {   free(cost);
            free(cost_prev);
            
            return sqrt(min_cost + cb[i+r+1]);
        }

        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return sqrt(final_dtw);
}




/** Unsigned arithmetic:
  * Given an 'index' and a 'window', get the start index corresponding to std::MAX(0, index-window) */
inline std::size_t cap_start_index_to_window(std::size_t index, std::size_t window){
    if(index>window){ return index-window; } else { return 0; }
}

/** Unsigned arithmetic:
  * Given an 'index', a 'window' and an 'end', get the stop index corresponding to std::MIN(end, index+window+1).
  * The expression index+window+1 is illegal for any index>0 as window could be MAX-1
  * */
inline std::size_t cap_stop_index_to_window_or_end(std::size_t index, std::size_t window, std::size_t end){
    // end-window is valid when window<end
    if(window<end && index+1<end-window){ return index + window + 1; } else { return end; }
}

double MON_dtw(
        const double* lines,
        const double* cols,
        const double* cb,
        int l,
        int w,
        double bsf
)
{
    // 1) --- Create the upper bound from bsf and a margin
    double UB= bsf - cb[w+1];

    // 2) --- Alias in line/column concept: we only allocate for the columns, using the smallest possible dimension.
    const size_t nbcols = l;
    const size_t nblines = l;

    // 3) --- Cap the windows.
    if (w > nblines) { w = nblines; }

    // 4) --- Buffers allocations
    // Add an extra column for the "matrix border" condition, init to +INF.
    // Using an unique contiguous array. Base indices are:
    // 'c' for the current row,
    // 'p' for the previous one
    std::vector<double> buffers_v((1+nbcols) * 2, POSITIVE_INFINITY);
    double *buffers = buffers_v.data();
    size_t c{1}, p{nbcols+2};                 // Account for the extra column (in front)

    // 5) --- Computation of DTW
    buffers[c-1] = 0;
    size_t next_start{0};
    size_t pruning_point{0};

    for(size_t i=0; i<nblines; ++i) {
        // --- --- --- --- Swap and variables init
        std::swap(c, p);
        if(i+w < nbcols){ UB = bsf - cb[i+w+1]; }
        const double li = lines[i];
        const std::size_t jStop = cap_stop_index_to_window_or_end(i, w, nbcols);
        //i + r + 1  m
        const std::size_t jStart = MAX(cap_start_index_to_window(i, w), next_start);
        //i-r 0
        std::size_t next_pruning_point = jStart; // Next pruning point init at the start of the line
        std::size_t j = jStart;
        next_start = jStart;

        // --- --- --- --- Init the first column
        buffers[c+j-1] = POSITIVE_INFINITY;
        double cost = POSITIVE_INFINITY;
        // --- --- --- --- Compute DTW up to the pruning point while advancing next_start: diag and top
        for(; j==next_start && j < pruning_point; ++j) {
            const auto d = DIST(li, cols[j]);
            cost = MIN(buffers[p + j - 1], buffers[p + j]) + d;
            buffers[c + j] = cost;
            if(cost<=UB){ next_pruning_point = j + 1;} else { ++next_start; }
        }

        // --- --- --- --- Compute DTW up to the pruning point without advancing next_start: prev, diag, top
        for(; j < pruning_point; ++j) {
            const auto d =  DIST(li, cols[j]);
            cost = MIN(cost, MIN(buffers[p + j - 1], buffers[p + j])) + d;
            buffers[c + j] = cost;
            if(cost<=UB){ next_pruning_point = j + 1;}
        }


        //
        // --- --- --- --- Compute DTW at "pruning_point": 2 cases
        if(j<jStop){
            const auto d = DIST(li, cols[j]);
            if(j==next_start){ // Advancing next start: only diag. Done if v>UB.
                cost = buffers[p + j - 1] + d;
                buffers[c + j] = cost;
                if(cost<=UB){ next_pruning_point = j + 1;}
                else
                {return POSITIVE_INFINITY; }
            } else { // Not advancing next start: at least a path possible in previous cells.
                cost =MIN(cost, buffers[p + j - 1]) + d;
                buffers[c + j] = cost;
                if(cost<=UB){ next_pruning_point = j + 1;}
            }
            ++j;
        } else if(j==next_start)
        { return POSITIVE_INFINITY; }


        // --- --- --- --- Compute DTW after "pruning_point": prev. Go on while we advance the next pruning point.
        for(;j==next_pruning_point && j<jStop;++j){
            const auto d =  DIST(li, cols[j]);
            cost = cost + d;
            buffers[c + j] = cost;
            if(cost<=UB){ ++next_pruning_point; }
        }
        // --- --- --- --- Row done, update the pruning point variable and the upper bound
        pruning_point=next_pruning_point;
    }// End for i loop


    // 6) --- If the pruning_point did not reach the number of columns, we pruned something
    if(pruning_point != nbcols)
    { return POSITIVE_INFINITY; }
    else {
        return sqrt(buffers[c+nbcols-1]);
    }
}
