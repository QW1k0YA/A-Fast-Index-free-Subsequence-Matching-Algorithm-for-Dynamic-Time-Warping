
#include "funcs.h"
#include "def.h"
#include "vector"
#include "fstream"
#include<iostream>
#include <iomanip>
#include "fftw3.h"
using namespace std;
double lb_sum = 0;
double lb_sum_ = 0;
double lb_kk_x_sum = 0;
double lb_kk_ul_sum = 0;
double lb_kk_called_num = 0;
double lb_kk_called_num_ = 0;
double lb_pe_num = 0;
double lb_t_called_num = 0;
double lb_q_called_num = 0;
double lb_q_fast_call_num=0;
/// 1. Given that lb_q has been thoroughly annotated, lb_q_fast, lb_t_new, and lb_data,
/// which are similar to lbq, will not be elaborated further here.
/// 2. There is a simple explanation of the functions related to Fourier Transform and Inverse Fourier Transform in sip.cpp.
/// 3. Please read the README.md before use code.
int main(int argc,char* argv[])
{

    double t1,t2;
    t1 = clock();
    vector<double> v_right;
    vector<double> v_left;
    double dtw_called_num = 0;
    double dtw_percent;
    double lb_q_n = 0;
    double lb_t_n = 0;
    double lb_q_n_new=0;
    double lb_d_n = 0;

    double threshold;
    vector<size_t> A;

    FILE *fp;            
    FILE *qp;            
    double bsf;          

    double d;
    long long i ,j;
    double ex_q , ex2_q,ex_t,ex2_t , mean, std;
    int m=-1, r,m_k = -1;

    double R;
    ///parameters
    m = atoi(argv[1]);
    R = atof(argv[2]);
    threshold = atof(argv[3]);
    int q_n = atoi(argv[4]);

    m_k = 2048;
    while(m_k <= 4*m)
    {
        m_k = m_k*2;
    }

    ///count the value of r
    if(R*m < 1)
    {
        r = 1;
    }
    else
    {
        r = round(R*m);
    }

    double threshold_2 = threshold*threshold;

    if(m<=8){
        cout<<"The length of query should be greater than 8!!!"<<endl;
        exit(-1);
    }

    long long k_m_l = m_k - m + 1;

    double *t = (double *) malloc(sizeof(double) * m_k);
    double *q = (double*) malloc(sizeof(double)*m);
    double *q_z = (double*) malloc(sizeof(double)*m);

    double *proj = (double *) malloc(sizeof(double) * m);
    double *u_p = (double*) malloc(sizeof(double)*m);
    double *l_p = (double*) malloc(sizeof(double)*m);

    bsf = INF;
    i = 0;
    j = 0;
    ex_q = ex2_q = 0;
    int count = 0;

    ///choose the query
    switch (q_n)
    {
        case 0:
        {
            qp = fopen(Q0,"r");
            break;
        }
        case 1:
        {
            qp = fopen(Q1,"r");
            break;
        }
        case 2:
        {
            qp = fopen(Q2,"r");
            break;
        }
        case 3:
        {
            qp = fopen(Q3,"r");
            break;
        }
        case 4:
        {
            qp = fopen(Q4,"r");
            break;
        }
        case 5:
        {
            qp = fopen(Q5,"r");
            break;
        }
        case 6:
        {
            qp = fopen(Q6,"r");
            break;
        }
        case 7:
        {
            qp = fopen(Q7,"r");
            break;
        }
        case 8:
        {
            qp = fopen(Q8,"r");
            break;
        }
        case 9:
        {
            qp = fopen(Q9,"r");
            break;
        }
        default:
        {
            exit(555);
            break;
        }
    }
    while(fscanf(qp,"%lf",&d) != EOF && i < m)
    {
        ex_q += d;
        ex2_q += d*d;
        q[i] = d;
        i++;
    }
    fclose(qp);

    mean = ex_q/m;
    std = ex2_q/m;
    std = sqrt(std-mean*mean);
    fflush(stdout);

    ///nomalise q
    for( i = 0 ; i < m ; i++ )
        q_z[i] = (q[i] - mean)/std;

    double q_max_cut = -1 * INF;
    double q_min_cut = INF;
    for(i = 0;i < m;i++)
    {
        if(q_z[i] < q_min_cut)
        {
            q_min_cut = q_z[i];
        }
        if(q_z[i] > q_max_cut)
        {
            q_max_cut = q_z[i];
        }
    }

    const double Q_min = q_min_cut;
    const double Q_max = q_max_cut;
    double derta = (Q_max - Q_min)/20;
    double inv_derta = 1/derta;
    i =0;
    int num = 1;
    double *lb_q = (double *) malloc(sizeof(double) * k_m_l);
    memset(lb_q,0,sizeof(double)*k_m_l);
    double *lb_t = (double *) malloc(sizeof(double) * k_m_l);
    double *lb_d = (double *) malloc(sizeof(double) * k_m_l);;

    int flag = 0;
    int flagofeof = 1;

    streampos sizeofT;
    
    string t_name = BI_T;
    ifstream input(t_name.c_str(),ios::binary | ios::in);
    input.seekg(0, ios::end);
    sizeofT = input.tellg();

    input.seekg(0, ios::beg);

    long long  ts_n=(sizeofT/(sizeof(double))) - m + 1;
    long long  forsize = (ts_n)/(k_m_l);
    int flag_of_tail = 0;
    size_t tail = (ts_n - m + 1)%(k_m_l);
    if(tail==0)
    {
        flag_of_tail = 1;
    }

    double *z  = (double *) malloc(sizeof(double) * m_k);

    fftw_complex  * XX = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)  * m_k);
    fftw_complex  * YY = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_k);
    fftw_complex  * YY_M = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_k);
    fftw_complex  * YY_UVM = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_k);
    fftw_complex  * YY_UVM_FAST = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_k);

    fftw_complex  * Z = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_k);
    fftw_complex  * ZZ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)  * m_k);
    fftw_complex  * q_PAD_FFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)  * m_k);
    fftw_complex  * QQ_PAD_FFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)  * m_k);

    fftw_plan p_forward = fftw_plan_dft_r2c_1d(m_k, z, ZZ,  FFTW_ESTIMATE);  
    fftw_plan p_backward = fftw_plan_dft_c2r_1d(m_k, ZZ, z, FFTW_ESTIMATE);

    double *QQ = (double *) malloc(sizeof(double) * m);
    ele_multi(q_z,q_z,QQ,m);
    
    double pad_value=0.1;

    inv_and_padding_y(m_k,m-6,QQ+3,QQ_PAD_FFT,p_forward);     
    inv_and_padding_y(m_k,m-6,q_z+3,q_PAD_FFT,p_forward);

    double *u_q = (double *) malloc(sizeof(double) * m);
    double *l_q = (double *) malloc(sizeof(double) * m);
    lower_upper_lemire(q_z,m,r,l_q,u_q);

    ///sort the order of q_z according to absolute values and envelopes respectively
    Index Q_tmp[m-6];
    Index Q_tmp_[m-6];
    for( i = 0; i < m-6; i++)
    {
        Q_tmp[i].value = q_z[i+3];
        Q_tmp[i].index = i+3;
    }
    for( i = 0; i < m-6; i++)
    {
        Q_tmp_[i].value = normal_cdf(u_q[i+3]+pad_value) - normal_cdf(l_q[i+3]-pad_value);
        Q_tmp_[i].index = i+3;
    }

    qsort(Q_tmp, m-6, sizeof(Index),comp);
    qsort(Q_tmp_, m-6, sizeof(Index),comp_);

    int *order_by_absQ =  (int *) malloc(sizeof(int) * m);
    double *q_zo = (double *) malloc(sizeof(double) * m);
    double *u_qo= (double *) malloc(sizeof(double) * m);
    double *l_qo= (double *) malloc(sizeof(double) * m);

    for( i=6; i<m; i++)
    {
        int o = Q_tmp[i-6].index;
        order_by_absQ[i] = o;
        q_zo[i] = q_z[o];
    }

    order_by_absQ[0] = 0;     order_by_absQ[1] = 1 ;     order_by_absQ[2] = 2;
    order_by_absQ[3] = m - 3;     order_by_absQ[4] = m - 2;    order_by_absQ[5] = m - 1;
    q_zo[0] = q_z[0];     q_zo[1] = q_z[1];     q_zo[2] = q_z[2];
    q_zo[3] = q_z[m-3];     q_zo[4] = q_z[m-2];     q_zo[5] = q_z[m-1];

    long long *order_by_envQ=  (long long *) malloc(sizeof(long long) * m);
    double *u_qo_= (double *) malloc(sizeof(double) * m);
    double *l_qo_= (double *) malloc(sizeof(double) * m);

    for( i=6; i<m; i++)
    {   int o = Q_tmp_[i-6].index;
        order_by_envQ[i] = o;
        u_qo_[i] = u_q[o];
        l_qo_[i] = l_q[o];
    }

    order_by_envQ[0] = 0;    order_by_envQ[1] = 1 ;     order_by_envQ[2] = 2;
    order_by_envQ[3] = m - 3;     order_by_envQ[4] = m - 2;     order_by_envQ[5] = m - 1;

    u_qo_[0] = u_q[0];     l_qo_[0] = l_q[0];
    u_qo_[1] = u_q[1];    l_qo_[1] = l_q[1];
    u_qo_[2] = u_q[2];    l_qo_[2] = l_q[2];
    u_qo_[3] = u_q[m -3];    l_qo_[3] = l_q[m -3];
    u_qo_[4] = u_q[m -2];    l_qo_[4] = l_q[m -2];
    u_qo_[5] = u_q[m -1];    l_qo_[5] = l_q[m -1];

    double M[m];
    
    for(int p = 0;p < m;p++)
    {
        if(normal_cdf(u_q[p]+pad_value) - normal_cdf(l_q[p]-pad_value) < 0.5) 
            M[p] = 1;
        else{
            M[p] = 0;
        }
    }

    int cnt_1 = 0;
    for(int ww = 0;ww < m;ww ++){
        cnt_1 += M[ww];
    }
    if(cnt_1 < 0.2*m)
    {
        int added_count=0;
        for(int p = 0;p < m;p++)
        {
            
            double local_pad_value = threshold/sqrt(0.2*m); 
            local_pad_value=min(local_pad_value,0.05);
            if(normal_cdf(u_q[p]+local_pad_value) - normal_cdf(l_q[p]-local_pad_value) < 0.5) 
            {
                M[p] = 1;
                added_count++;
            }
        }
        cout<<"less than 0.1m of 1, added : "<< added_count<< "      p: "<< cnt_1*1.0/m  <<endl;
    }
   
    cnt_1 = 0;
    for(int ww = 0;ww < m;ww ++)  cnt_1 += M[ww];

    M[0]=0; M[1]=0; M[2]=0; M[m-1]=0; M[m-2]=0; M[m-3]=0;

    double* UM_lbq = NULL;
    double* VM_lbq = NULL;
    double* UU = NULL;
    double* VV = NULL;
    double* UV_PLUS_M = NULL;

    UM_lbq = (double*)malloc(m * sizeof(double));
    VM_lbq = (double*)malloc(m * sizeof(double));
    UU = (double*)malloc(m * sizeof(double));
    VV = (double*)malloc(m * sizeof(double));
    UV_PLUS_M = (double*)malloc(m * sizeof(double));

    ele_multi(M,l_q,VM_lbq,m);
    ele_multi(M,u_q,UM_lbq,m);
    ele_plus(UM_lbq,VM_lbq,UV_PLUS_M,m);

    double M_FAST[m];
    for(int mm = 0; mm < m;mm++)
    {
        M_FAST[mm] = 1;
    }
    M_FAST[0]=0; M_FAST[1]=0; M_FAST[2]=0; M_FAST[m-1]=0; M_FAST[m-2]=0; M_FAST[m-3]=0;

    double *UM_lbq_FAST = (double *) malloc(sizeof(double) * m);
    double *VM_lbq_FAST = (double *) malloc(sizeof(double) * m);
    double *UV_PLUS_M_FAST = (double *) malloc(sizeof(double) * m);

    ele_multi(M_FAST,l_q,VM_lbq_FAST,m);
    ele_multi(M_FAST,u_q,UM_lbq_FAST,m);
    ele_plus(UM_lbq_FAST,VM_lbq_FAST,UV_PLUS_M_FAST,m);

    ele_multi(u_q,u_q,UU,m);
    ele_multi(l_q,l_q,VV,m);

    inv_and_padding_y(m_k,m,M,YY_M,p_forward);
    inv_and_padding_y(m_k,m,UV_PLUS_M,YY_UVM,p_forward);
    inv_and_padding_y(m_k,m,UV_PLUS_M_FAST,YY_UVM_FAST,p_forward);

    double* u_t = (double *) malloc(sizeof(double) * m);
    double* l_t = (double *) malloc(sizeof(double) * m);
    double* u_t_test = (double *) malloc(sizeof(double) * m);
    double* l_t_test = (double *) malloc(sizeof(double) * m);
    double* u_t_ = (double *) malloc(sizeof(double) * m_k);
    double* l_t_ = (double *) malloc(sizeof(double) * m_k);

    int i_ = 0;
    double lbq_target = ((m_k-m+1)*(1-log(m_k)/m));

    double lb_target2= (m_k-m+1)*(0.8);
    bool  lbf_target = false;

    double  M_module_2 = 0;
    double  M_module_2_FAST = 0;
    for(i = 0;i < m;i++)
    {
        M_module_2 += M[i];
        M_module_2_FAST += M_FAST[i];
    }

    if(M_module_2/m <=  0.8)
    {
        lbf_target = true;
    }

    double c1 = dot_mul(UU,M,m);
    double c2 = dot_mul(VV,M,m);
    double c3 = dot_mul(u_q,M,m);
    double c4 = dot_mul(l_q,M,m);
    double c12 = c1 + c2;
    double c34 = c3 + c4;

    double c1_fast = dot_mul(UU,M_FAST,m);
    double c2_fast = dot_mul(VV,M_FAST,m);
    double c3_fast = dot_mul(u_q,M_FAST,m);
    double c4_fast = dot_mul(l_q,M_FAST,m);
    double c12_fast = c1_fast + c2_fast;
    double c34_fast = c3_fast + c4_fast;

    double* share_k_m_l1 = (double *) malloc(sizeof(double) * (m_k)); 
    double* share_k_m_l2 = (double *) malloc(sizeof(double) * (m_k)); 
    double* share_k_m_l3 = (double *) malloc(sizeof(double) * (m_k)); 
    double* share_k_m_l4 = (double *) malloc(sizeof(double) * (m_k)); 

    double *share_m_k1 = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k2 = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k3 = (double *) malloc(sizeof(double) * m_k);

    double *share_m_k_1 = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k_2 = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k_3 = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k_4_TTonly = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k_5 = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k_6 = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k_7 = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k_8 = (double *) malloc(sizeof(double) * m_k);
    double *share_m_k_9 = (double *) malloc(sizeof(double) * m_k);

    double *cb = (double *) malloc(sizeof(double) * m);;
    memset(cb, 0, sizeof(double) * m);
    double* cb1 = (double *) malloc(sizeof(double) * m);;
    memset(cb1, 0, sizeof(double) * m);
    double *cb2 = (double *) malloc(sizeof(double) * m);;
    memset(cb2, 0, sizeof(double) * m);

    double **shared_mk_2_vector = (double **)malloc(m_k * sizeof(double *));
    if (shared_mk_2_vector == NULL) {
        fprintf(stderr, "fail in malloc\n");
        return 1;
    }

    for (i = 0; i < m_k; i++) {
        shared_mk_2_vector[i] = (double *)malloc(2 * sizeof(double));
        if (shared_mk_2_vector[i] == NULL) {
            fprintf(stderr, "fail in malloc\n");
            
            for (j = 0; j < i; j++) {
                free(shared_mk_2_vector[j]);
            }
            free(shared_mk_2_vector);
            return 1;
        }
    }

    double *miu_t = (double*) malloc(sizeof(double)*k_m_l);
    double *si_t = (double*) malloc(sizeof(double)*k_m_l);
    double *tab_q = (double*) malloc(sizeof(double)*m*20);
    table_of_q(q_z, m, r, tab_q, Q_min, derta);

    size_t bin_num = BIN;
    int bin[bin_num];
    memset(bin, 0, bin_num * sizeof(int));
    q_max_cut=min(Q_max, 5.0); q_min_cut=max(Q_min, -5.0);
    double bin_num_divide_Q_max_Q_min = (bin_num) / (q_max_cut - q_min_cut);
    for(i = 0;i < m;i++)
    {
        bin[H(q_z[i], bin_num_divide_Q_max_Q_min, q_min_cut, bin_num)] += 1;
    }
    for(i = 1;i < bin_num;i ++)
    {
        bin[i] += bin[i-1];
    }

    double time_for_now=clock();
    cout << "the prepare takes : " << (time_for_now - t1)/CLOCKS_PER_SEC << "sec"<< endl;

    double lbq_cnt = 0;

    int forsize_divide_100 = forsize/100;

    double* M_data = (double *) malloc(sizeof(double) * m);
    double* MQ = (double *) malloc(sizeof(double) * m);
    double D_uv = ED(UM_lbq, VM_lbq, m);

    for(long long ll = 0; ll < forsize;ll++)
    {

        if(forsize_divide_100)
        {
            if((ll%(forsize_divide_100)) == 0)
            {
                cout << "\b\b\b\b";

                cout << "|" << setw(3) << i_ << "%";

                i_++;
            }
        }

        ///If it is not the first time reading the data, then read a portion of the overlapping section.
        if(ll!= 0)
        {
            memcpy(t,t+k_m_l,sizeof(double)*m);
        }
        else
        {
            input.read((char *)t, (m-1)* sizeof(double));
        }

        input.read((char *)(t + m -1), (k_m_l)* sizeof(double));

        int lb_n = 0, lb_n_new=0;
        bool flag_lbq = 0;
        mvmean(t,m_k,m,miu_t,si_t);

        lb_n = LB_q(q_z, t, m_k, m, lb_q, threshold, miu_t, si_t, XX, YY, ZZ, p_forward, p_backward, YY_M, YY_UVM,
                    share_m_k1, share_m_k2, share_m_k3, M_module_2, c12, c34,
                    shared_mk_2_vector, share_m_k_4_TTonly, D_uv);
        lb_q_called_num ++;
        lb_q_n+=lb_n;

        ///If lbq does not filter out enough sequences, then invoke lbqfast.
        if (lbf_target && (lb_n < lbq_target)){

            lb_n_new= LB_q_fast(q_z, t, m_k, m, lb_q, threshold, miu_t, si_t, Z, XX, ZZ, p_backward,
                            UM_lbq_FAST, VM_lbq_FAST, YY_UVM_FAST, share_m_k1, share_m_k2, share_m_k3, M_module_2_FAST, c12_fast, c34_fast,
                            shared_mk_2_vector, share_m_k_4_TTonly);
            
            lb_q_n_new+=lb_n_new;
            lb_q_fast_call_num++;
        }
        int flag_lemire = 0;


        ///Invoke lbdata and lbtnew for further filtering.
        if(lb_n+lb_n_new > lbq_target)
        {
            flag_lbq = 1;
        }
        else
        {
            lower_upper_lemire(t, m_k, r, l_t_, u_t_);
            int lb_n2= LB_t_new(m_k, m, lb_t, u_t_, l_t_, threshold, miu_t, si_t, XX, q_PAD_FFT, QQ_PAD_FFT, ZZ,
                                p_forward, p_backward,
                                share_m_k1, share_m_k2, share_m_k3, share_k_m_l1, share_k_m_l2, share_k_m_l3,
                                share_k_m_l4, share_m_k_1, share_m_k_2, share_m_k_3, share_m_k_4_TTonly, share_m_k_6,
                                share_m_k_9, Q_min, shared_mk_2_vector, lb_q, bin, bin_num_divide_Q_max_Q_min);
           
            lb_t_called_num ++;
            lb_t_n+=lb_n2;
            flag_lemire = 0;

            if ( lb_n+lb_n_new+lb_n2 <= lb_target2){
                int mlen;
                
                 mlen = 0.1 * m;
                MASK_DATA(order_by_absQ,m,M_data,mlen);

                lb_d_n += lb_data(lb_d, M_data, q_z, m, u_t_, l_t_, miu_t, MQ, m_k, XX, p_forward, p_backward, share_m_k1,
                                  share_m_k2,
                                  share_m_k3, share_m_k_1, YY, ZZ, share_m_k_2, share_k_m_l1, share_k_m_l2, si_t,
                                  share_m_k_3,
                                  share_m_k_5, share_m_k_6, share_m_k_7, share_m_k_8, threshold_2, shared_mk_2_vector,
                                  lb_q);
            }
       }

        double  lb_kk_x, lb_kk_uv;
        double u_t_ul[m];  double l_t_ul[m];
        
        for(int p = 0;p < k_m_l;p++)
        {
            double *t_temp = t + p;
            lb_sum += (lb_q[p]+1)*lb_q[p];

            if(threshold >= lb_q[p])
            {

                double lb_k = 0;
                double lb_k2 = 0;
                lb_kk_called_num ++;
                if(LB_KK_norm_X_on_the_fly_enhance(t + p, u_qo_, l_qo_, m, miu_t[p], si_t[p],
                                                   threshold_2, order_by_envQ, cb1, lb_k, tab_q,
                                                   shared_mk_2_vector[p], Q_min, inv_derta))
                {
                    if(flag_lbq)
                    {
                        ///count the envelope of t
                        if(flag_lemire > p + 2*r && flag_lemire)
                        {
                            lower_upper_lemire(t+flag_lemire-2*r, m - (flag_lemire - p) + 2*r, r, l_t_ul, u_t_ul);
                            memcpy(u_t_ + flag_lemire - r,u_t_ul + r,sizeof(double)*(m - (flag_lemire - p) + r));
                            memcpy(l_t_ + flag_lemire - r,l_t_ul + r,sizeof(double)*(m - (flag_lemire - p) + r));
                        }
                        else
                        {
                            lower_upper_lemire(t+p, m, r, l_t_+p, u_t_+p);
                            flag_lemire = p + m;
                        }
                    }
                    lb_kk_called_num_ ++;

                    if(LB_KK_norm_UL_on_the_fly(q_zo, u_t_+p, l_t_+p, m, miu_t[p], si_t[p],
                                                threshold_2, order_by_absQ, cb2,lb_k2,shared_mk_2_vector[p]))
                    {

                        double t_temp_[m];
                        for(int tt = 0;tt < m;tt++)
                        {
                            t_temp_[tt] = (t_temp[tt] - miu_t[p]) / si_t[p];
                        }
                        
                        lb_pe_num++;
                        ///lbPetitjean_new provide a tighter lower bound
                        if(lbPetitjean_new(proj, q_z, u_p, l_p, u_q, l_q, t_temp_, u_t_ + p, l_t_ + p, r, m,
                                           threshold_2, lb_k2,
                                           cb2, miu_t[p], si_t[p]))

                        {
                            dtw_called_num++;
                            if (lb_k > lb_k2)
                            {
                                cb[m-1]=cb1[m-1];
                                for(int k=m-2; k>=0; k--)
                                    cb[k] = cb[k+1]+cb1[k];
                            }
                            else
                            {
                                cb[m-1]=cb2[m-1];
                                for(int k=m-2; k>=0; k--)
                                    cb[k] = cb[k+1]+cb2[k];
                            }

                            if(threshold >= dtw(q_z, t_temp_, m, r, threshold_2, cb))
                            {
                                
                                A.push_back(ll*k_m_l  + p);
                                
                            }
                        }
                    }

                }
            }

        }
    }

    ///Handle the remaining part
    if(!flag_of_tail)
    {
        for(int f=0; f<m - 1; f++)
            t[f] = t[k_m_l+f];
        input.read((char *)(t + m -1), (tail)* sizeof(double));

        int siz_t = tail + m -1;
        size_t st = sizeofT/(sizeof(double)) - siz_t;

        mvmean(t, siz_t, m, miu_t, si_t);
        lower_upper_lemire(t, siz_t, r, l_t_, u_t_ );

        for(int p = 0;p < tail;p++)
        {
            double *t_temp = t + p;
            double lb_Kim=  LB_KIM(t, q_z, miu_t, si_t, threshold, m, shared_mk_2_vector[p]);
            if (lb_Kim> threshold_2)      continue;

            double lb_k = 0;
            double lb_k2 = 0;
            lb_kk_called_num++;
            if(si_t[p] == 0 || isnan(si_t[p]))
                continue;

            if(LB_KK_norm_X_on_the_fly_enhance(t + p, u_qo_, l_qo_, m, miu_t[p], si_t[p], threshold_2, order_by_envQ,
                                               cb1, lb_k, tab_q, shared_mk_2_vector[p], Q_min, inv_derta))
            {
                lb_kk_called_num_++;
                if(LB_KK_norm_UL_on_the_fly(q_zo, u_t_ + p, l_t_ + p, m, miu_t[p], si_t[p], threshold_2, order_by_absQ, cb2, lb_k2, shared_mk_2_vector[p]))
                {
                    double t_temp_[m];
                    for(int tt = 0;tt < m;tt++)
                    {
                        t_temp_[tt] = (t_temp[tt] - miu_t[p]) / si_t[p];
                    }
                   
                    lb_pe_num++;
                    
                    {
                        dtw_called_num++;
                        if (lb_k > lb_k2)
                        {
                            cb[m-1]=cb1[m-1];
                            for(int k=m-2; k>=0; k--)
                                cb[k] = cb[k+1]+cb1[k];
                        }
                        else
                        {
                            cb[m-1]=cb2[m-1];
                            for(int k=m-2; k>=0; k--)
                                cb[k] = cb[k+1]+cb2[k];
                        }

                        if(threshold >= dtw(q_z, t_temp_, m, r, threshold_2, cb))
                        {
                            
                            A.push_back(st  + p);

                        }
                    }

                }
            }
        }
    }

    ///free all
    free(QQ);
    free(order_by_absQ);
    free(order_by_envQ);
    free(q_zo);
    free(u_qo_);
    free(l_qo_);
    free(UM_lbq);
    free(VM_lbq);
    free(UU);
    free(VV);
    free(UV_PLUS_M);
    free(z);
    free(miu_t);
    free(si_t);
    free(tab_q);

    free(lb_q);
    free(lb_t);
    free(u_q);

    free(l_q);
    free(t);
    free(q);
    free(q_z);

    free(cb);
    free(cb1);
    free(cb2);
    free(l_t_);
    free(u_t_);

    free(UM_lbq_FAST);
    free(VM_lbq_FAST);
    free(UV_PLUS_M_FAST);
    fftw_free(XX); fftw_free(YY);
    fftw_free(Z); fftw_free(ZZ);
    fftw_free(q_PAD_FFT);
    fftw_free(QQ_PAD_FFT);

    fftw_destroy_plan(p_backward);
    fftw_destroy_plan(p_forward);

    free(share_m_k1);
    free(share_m_k3);
    free(share_m_k2);

    free(share_k_m_l1);
    free(share_k_m_l2);
    free(share_k_m_l3);
    free(share_k_m_l4);

    free(share_m_k_1);
    free(share_m_k_2);
    free(share_m_k_3);
    free(share_m_k_4_TTonly);
    free(share_m_k_5);
    free(share_m_k_6);
    free(share_m_k_7);
    free(share_m_k_8);
    free(share_m_k_9);
    free(u_p);
    free(l_p);
    free(proj);
    free(u_t);
    free(l_t);

    input.close();
    t2 = clock();
    cout << "lbq_cnt = " << lbq_cnt << endl;
    cout << "total ex time is " << (t2 - t1)/CLOCKS_PER_SEC << "sec"<< endl;

    int A_size = A.size();
    cout << "A_size = " << A.size() <<endl;
    cout <<"input: " <<  m << " " << atof(argv[2]) << " " << atof(argv[3]) << " "<< atoi(argv[4]) << endl;
    cout << "lb_sum = " << lb_sum << endl;
    cout << "lb_sum_= " << lb_sum_ << endl;
    cout << "lb_kk_x_sum = " << lb_kk_x_sum << endl;
    cout << "lb_kk_uv_sum = " << lb_kk_ul_sum << endl;
    cout << "lb_q_called = " <<lb_q_called_num <<endl;
    cout <<"lb_q_f_called = "<<lb_q_fast_call_num<<endl;
    cout << "lb_t_called_num = " <<lb_t_called_num << endl;
    cout << "skipped lbt computation = " << 1- (lb_q_called_num)/(1.0*forsize) <<endl;
    
    cout <<"#################STATIC#####################"<<endl;
    cout << "lbq purns = " << lb_q_n/ts_n <<endl; 
    cout <<"lbq_fast purns = "<< lb_q_n_new/ts_n<<endl;
    cout << "lbt purns = " <<  lb_t_n/ts_n <<endl;
    cout << "lb_data purns = " <<  lb_d_n/ts_n <<endl;
    cout << "LBKKX purns = " << (lb_kk_called_num-lb_kk_called_num_)/(sizeofT/(sizeof(double)) - m + 1) <<endl;
    cout << "LBKKUL purns = " << (lb_kk_called_num_-lb_pe_num)/(sizeofT/(sizeof(double)) - m + 1) <<endl;
    cout << "LB_PETITJEAN = " << (lb_pe_num - dtw_called_num)/ (sizeofT/(sizeof(double)) - m + 1) <<endl;;
    cout << "dtw_percent = " << dtw_called_num/ts_n << endl;
    cout << endl;

}
