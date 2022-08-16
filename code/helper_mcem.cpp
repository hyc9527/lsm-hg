// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

/*
    MCEM helpers in cpp, including:
        (i)  estep
        (ii) likelihood functions
*/

static double const log2pi = std::log(2.0 * M_PI);

// Function 0.1 Obsolete
// Remark: since the function takes column vector only and thus arma::norm is prefered
// [[Rcpp::export]]
double norm_l2(arma::vec x)
{
    return arma::norm(x, 2);
}

// Function 0.2
// [[Rcpp::export]]
double size_triangle(arma::irowvec index, arma::mat Feat)
{
    double a, b, c, s;
    a = arma::norm(Feat.row(index[0]) - Feat.row(index[1]), 2); //Rcpp slicing index starts at 0
    b = arma::norm(Feat.row(index[0]) - Feat.row(index[2]), 2);
    c = arma::norm(Feat.row(index[1]) - Feat.row(index[2]), 2);
    if (abs(a) < 1e-10 || abs(b) < 1e-10 || abs(c) < 1e-10)
    {
        return 0;
    }
    else
    {
        s = 0.5 * (a + b + c);
    }
    return sqrt(s * (s - a) * (s - b) * (s - c));
}

// Function 0.3
//   Cited from  https://gallery.rcpp.org/articles/simulate-multivariate-normal/
//   each row of returned matrix is distributed as mvrnorm(mu, sigma)
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat Sigma)
{
    int ncols = Sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
}

/* C++ version of the dtrmv BLAS function */
// helper in dmvnrm_arma_fast()
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat)
{
    arma::uword const n = trimat.n_cols;

    for (unsigned j = n; j-- > 0;)
    {
        double tmp(0.);
        for (unsigned i = 0; i <= j; ++i)
            tmp += trimat.at(i, j) * x[i];
        x[j] = tmp;
    }
}

// https://gallery.rcpp.org/articles/dmvnorm_arma/
// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &sigma,
                           bool const logd = false)
{
    using arma::uword;
    uword const n = x.n_rows,
                xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())),
                 constants = -(double)xdim / 2.0 * log2pi,
                 other_terms = rootisum + constants;
    arma::rowvec z;
    for (uword i = 0; i < n; i++)
    {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);
        out(i) = other_terms - 0.5 * arma::dot(z, z);
    }
    if (logd)
        return out;
    return exp(out);
}

// Function 1.1
// Return volume of k-cvxhull ONLY for k = 1,2,3 !
// [[Rcpp::export]]
double exact_volume_cvxhull(arma::irowvec subset, arma::mat Feat)
{
    int k = subset.n_elem;
    double res;
    if (k == 1)
    {
        res = 0;
    }
    else if (k == 2)
    {
        res = arma::norm(Feat.row(subset[0]) - Feat.row(subset[1]), 2);
    }
    else
    {
        res = size_triangle(subset, Feat);
    }
    return res;
}

// Function 1.2' Call a r function from cpp.
//      Return volume of k-cvxhull ONLY for k \geq 4 !
//      https://stackoverflow.com/questions/38016851/call-r-functions-in-rcpp
//      Call r Function 1.2 'volume.cvxhull' from rcpp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! A SLOW WORKAROUND

// [[Rcpp::export]]
double volume_cvxhull_4pt_cpp_call(const arma::irowvec subset, const arma::mat Feat)
{

    // Obtain environment containing function
    //Rcpp::Environment func_env("R_GlobalEnv");
    // Make function callable from C++
    Rcpp::Function func("volume.cvxhull.4pt.plus"); // R fucntion 'volume.cvxhull.4pt.plus' is imported into G_GlobalEnv;
    Rcpp::NumericVector res = func(subset, Feat);   // example of additional param
    return res[0];
}

// Function 1.2
//      For k = 2,3 exact volume is returned
//      For k \geq 4, call a r function 'volume.cvxhull.4pt.plus' from cpp. COULD BE SLOW !!!
//      [[Rcpp::export]]
double find_volume_cvxhull(const arma::irowvec subset, const arma::mat Feat)
{
    int k = subset.n_elem;
    double vol;
    if (k <= 3)
    {
        vol = exact_volume_cvxhull(subset, Feat);
    }
    else
    {
        vol = volume_cvxhull_4pt_cpp_call(subset, Feat);
    }
    return vol;
}

// Function 3.1: a special case of r function 3.1
// [[Rcpp::export]]
double log_shape_den_gamma(double t, arma::vec param)
{
    return param[0] * log(param[1]) - lgamma(param[0]) + (param[0] - 1) * log(t) - param[1] * t;
}

//Function 3.2
// [[Rcpp::export]]
double log_likelihood_shape_cvxhull_gamma(const arma::ivec count, const arma::mat Feat, const arma::imat Subset, arma::vec param)
{
    int n_k = Subset.n_rows;
    arma::vec kern(n_k);
    double sum_kern = 0;
    for (int j = 0; j < n_k; j++) // calculate normalizing term
    {
        //Rf_PrintValue(Rcpp::wrap(Subset.row(j) - 1));
        double t_j = find_volume_cvxhull(Subset.row(j) - 1, Feat); // arg 'Subset' is a mat of all k-sets of [n], where elements are using R indexing style from 1:n
        if (t_j < 1e-300)
        {
            t_j = 1e-300; // dealing with almose linear dependent cases:
                          // set volume  to be 1e-100 if it is too small.
        }                 //    if volume is zero, log_shape_den_gamma() returns NAN
        //Rf_PrintValue(Rcpp::wrap(t_j));
        kern[j] = exp(log_shape_den_gamma(t_j, param));
        sum_kern += kern[j];
        //if (count[j] > 0)
        //#{
        //llkhd += count[j] * log(kern[j]);
        //#}
    }
    //Rf_PrintValue(Rcpp::wrap(kern.t()));
    //Rf_PrintValue(Rcpp::wrap(sum_kern));
    //int num_zero_count = 0;
    double llkhd = 0;
    int num_zeros_count = 0;
    for (int j = 0; j < n_k; j++)
    {
        if (count[j] == 0)
        {
            //num_zero_count += 1;
            llkhd += log(sum_kern - kern[j]);
            num_zeros_count += 1;
            //Rf_PrintValue(Rcpp::wrap(log(sum_kern - kern[j])));
        }
        else
        {
            llkhd += count[j] * log(kern[j]);
            //Rf_PrintValue(Rcpp::wrap(log(kern[j])));
        }
    }
    llkhd -= (sum(count) + num_zeros_count) * log(sum_kern);
    ////Rf_PrintValue(llkhd);
    //Rf_PrintValue(llkhd);
    return llkhd;
}

// Function 3.12
//       R verions is adopted from https: //rpubs.com/mengxu/procrustes_analysis
//       return a matrix(i) of which cols are of zero mean;
//       (ii).|| mat || _f = n * p
// [[Rcpp::export]]
arma::mat standarize_matrix(const arma::mat X)
{
    int n = arma::size(X)[0];
    int p = arma::size(X)[1];
    arma::mat res(n, p);
    for (int j = 0; j < p; j++)
    {
        double col_mean = arma::mean(X.col(j));
        res.col(j) = X.col(j) - col_mean;
    }
    double norm_f = norm(res, "fro");
    res = res / norm_f * n * p;
    return res;
}

// Function 3.13
//        R version, adopted from https://rpubs.com/mengxu/procrustes_analysis
// [[Rcpp::export]]
arma::mat find_procrustean_transform(const arma::mat std_Z0, const arma::mat std_Z)
{
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::mat X = std_Z * std_Z0.t();
    arma::svd(U, s, V, X);
    arma::mat res = V * U.t() * std_Z;
    return res;
}

// Function:  original estep  with a bunch of undeadly bugs which unluckily works
// directly translated from a R version mcem.estep, is labeled as a<1>b<1>c<1>
// [[Rcpp::export]]
Rcpp::List mcem_estep_cpp(const int em_round, const int min_sweeps,
                          arma::mat Sigma_prop,
                          Rcpp::List res_from_last_step,
                          const int tot_author, const arma::imat all_subset_mat, const arma::irowvec data_s,
                          bool if_message_reported)
{
    int rej = 0;
    Rcpp::List beta_lst = res_from_last_step["beta.lst"];
    int len_beta_lst = beta_lst.size();
    double beta_last = beta_lst[len_beta_lst - 1];
    arma::vec param(2, arma::fill::ones);
    param[1] = beta_last;
    arma::mat X = res_from_last_step["last.feat.mat"];
    arma::mat X_pivot = X;
    int num_sweeps = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List llkhd_mstep = res_from_last_step["likelihood.mstep.lst"];
    arma::vec old_llkhd(1);
    old_llkhd[0] = llkhd_mstep[em_round - 2];
    int m_t = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List Feat_sweep_lst(m_t);
    //Rprintf("Estep is starting.\n");
    for (arma::uword j = 0; j < m_t; j++)
    {
        for (int i = 0; i < tot_author; i++)
        {
            //Rprintf("Iter (%i, %i) \n", j, i);
            arma::rowvec x_old = X.row(i);
            arma::mat x_new = mvrnormArma(1, x_old.t(), Sigma_prop);
            x_new = arma::conv_to<arma::rowvec>::from(x_new);
            arma::mat X_new = X;
            X_new.row(i) = x_new;
            X_new = standarize_matrix(X_new);
            x_new = X_new.row(i);
            arma::vec curr_llkhd(1);
            curr_llkhd = log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                            X_new,
                                                            all_subset_mat,
                                                            param);
            arma::vec log_ratio(1);
            log_ratio = curr_llkhd;
            arma::mat x_new_in_mat = arma::conv_to<arma::mat>::from(x_new);
            arma::mat x_old_in_mat = arma::conv_to<arma::mat>::from(x_old);
            log_ratio += dmvnrm_arma_fast(x_new_in_mat, x_old, Sigma_prop, true); // prior density in numerator  // BUG!!!
            log_ratio -= old_llkhd;
            log_ratio -= dmvnrm_arma_fast(x_old_in_mat, x_old, Sigma_prop, true); //  prior density in denominator // BUG!!!
            bool valid = all(log(arma::vec(1, arma::fill::randu)) < log_ratio);
            if (valid)
            {
                X = X_new;
                old_llkhd = curr_llkhd;
            }
            else
            {
                rej += 1;
            }
        }
        //Rprintf("Sweep no. %i is done.\n", j);
        //Rcpp::print(Rcpp::wrap(X_pivot));
        //Rcpp::print(Rcpp::wrap(X));
        arma::mat X_prcr = find_procrustean_transform(X_pivot, X);
        //Rcpp::print(Rcpp::wrap(X_prcr));
        Feat_sweep_lst[j] = X_prcr;
    }
    int warm_start = ceil(m_t / 2);
    Rcpp::List Feat_sweep_lst_warm = Feat_sweep_lst[Rcpp::Range(warm_start, m_t - 1)];
    arma::vec likelihood_estep(1, arma::fill::zeros);
    for (int t = 0; t < Feat_sweep_lst.size(); t++)
    {
        likelihood_estep += log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                               Feat_sweep_lst[t],
                                                               all_subset_mat,
                                                               param);
    }
    likelihood_estep = likelihood_estep / Feat_sweep_lst.size();
    Rcpp::List res = Rcpp::List::create();
    res["Feat.lst"] = Feat_sweep_lst;
    res["likelihood.estep"] = likelihood_estep;
    res["last.feat.mat"] = X;
    res["rej"] = rej;
    return (res);
}

// Function estep with variation
// [[Rcpp::export]]
Rcpp::List mcem_estep_cpp_a1b2c1(const int em_round, const int min_sweeps,
                                 arma::mat Sigma_prop,
                                 Rcpp::List res_from_last_step,
                                 const int tot_author, const arma::imat all_subset_mat, const arma::irowvec data_s,
                                 bool if_message_reported)
{
    int rej = 0;
    Rcpp::List beta_lst = res_from_last_step["beta.lst"];
    int len_beta_lst = beta_lst.size();
    double beta_last = beta_lst[len_beta_lst - 1];
    arma::vec param(2, arma::fill::ones);
    param[1] = beta_last;
    arma::mat X = res_from_last_step["last.feat.mat"];
    arma::mat X_pivot = X;
    int num_sweeps = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List llkhd_mstep = res_from_last_step["likelihood.mstep.lst"];
    arma::vec old_llkhd(1);
    old_llkhd[0] = llkhd_mstep[em_round - 2];
    arma::rowvec mu_0 = arma::rowvec(2, arma::fill::zeros);
    arma::mat Sigma_0 = arma::mat(2, 2, arma::fill::eye);
    int m_t = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List Feat_sweep_lst(m_t);
    for (arma::uword j = 0; j < m_t; j++)
    {
        for (int i = 0; i < tot_author; i++)
        {
            arma::rowvec x_old = X.row(i);
            arma::mat x_new = mvrnormArma(1, x_old.t(), Sigma_prop);
            x_new = arma::conv_to<arma::rowvec>::from(x_new);
            arma::mat X_new = X;
            X_new.row(i) = x_new;
            //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& variation a<1>
            X_new = standarize_matrix(X_new);
            x_new = X_new.row(i);
            //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            arma::vec curr_llkhd(1);
            curr_llkhd = log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                            X_new,
                                                            all_subset_mat,
                                                            param);
            arma::vec log_ratio(1);
            log_ratio = curr_llkhd;
            arma::mat x_new_in_mat = arma::conv_to<arma::mat>::from(x_new);
            arma::mat x_old_in_mat = arma::conv_to<arma::mat>::from(x_old);
            //log_ratio += dmvnrm_arma_fast(x_new_in_mat, mu_0, Sigma_prop, true);  // BUG !!!
            //log_ratio += dmvnrm_arma_fast(x_new_in_mat, x_old, Sigma_prop, true) ;   // proposal density for x_new //BUG!!!
            log_ratio += dmvnrm_arma_fast(x_new_in_mat, mu_0, Sigma_0, true); // prior density in numerator  // CHEATING !!!
            log_ratio -= old_llkhd;
            //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& variation b<2>
            log_ratio -= dmvnrm_arma_fast(x_old_in_mat, mu_0, Sigma_0, true); // prior density for x_old in denominator // CHEATING !!!
            //Rcpp::print(Rcpp::wrap(dmvnrm_arma_fast(x_new_in_mat, mu_0, Sigma_0, true))); // test
            //Rcpp::print(Rcpp::wrap(dmvnrm_arma_fast(x_old_in_mat, mu_0, Sigma_0 , true)));    // test
            //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            bool valid = all(log(arma::vec(1, arma::fill::randu)) < log_ratio);
            if (valid)
            {
                X = X_new;
                old_llkhd = curr_llkhd;
            }
            else
            {
                rej += 1;
            }
        }
        //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& variation c<1>
        arma::mat X_prcr = find_procrustean_transform(X_pivot, X);
        Feat_sweep_lst[j] = X_prcr;
        //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    }
    int warm_start = ceil(m_t / 2);
    Rcpp::List Feat_sweep_lst_warm = Feat_sweep_lst[Rcpp::Range(warm_start, m_t - 1)];
    arma::vec likelihood_estep(1, arma::fill::zeros);
    for (int t = 0; t < Feat_sweep_lst.size(); t++)
    {
        likelihood_estep += log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                               Feat_sweep_lst[t],
                                                               all_subset_mat,
                                                               param);
    }
    likelihood_estep = likelihood_estep / Feat_sweep_lst.size();
    Rcpp::List res = Rcpp::List::create();
    res["Feat.lst"] = Feat_sweep_lst;
    res["likelihood.estep"] = likelihood_estep;
    res["last.feat.mat"] = X;
    res["rej"] = rej;
    return (res);
}

// Function estep_2 in mcem
// [[Rcpp::export]]
Rcpp::List mcem_estep_cpp_a1b2c3(const int em_round, const int min_sweeps,
                                 arma::mat Sigma_prop,
                                 Rcpp::List res_from_last_step,
                                 const int tot_author, const arma::imat all_subset_mat, const arma::irowvec data_s,
                                 bool if_message_reported)
{
    int rej = 0;
    Rcpp::List beta_lst = res_from_last_step["beta.lst"];
    int len_beta_lst = beta_lst.size();
    double beta_last = beta_lst[len_beta_lst - 1];
    arma::vec param(2, arma::fill::ones);
    param[1] = beta_last;
    arma::mat X = res_from_last_step["last.feat.mat"];
    arma::mat X_pivot = X;
    int num_sweeps = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List llkhd_mstep = res_from_last_step["likelihood.mstep.lst"];
    arma::vec old_llkhd(1);
    old_llkhd[0] = llkhd_mstep[em_round - 2];
    int m_t = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List Feat_sweep_lst(m_t);
    arma::mat X_anchor = res_from_last_step["anchor.feat.mat"];
    arma::rowvec mu_0 = arma::rowvec(2, arma::fill::zeros);
    arma::mat Sigma_0 = arma::mat(2, 2, arma::fill::eye);
    for (arma::uword j = 0; j < m_t; j++)
    {
        for (int i = 0; i < tot_author; i++)
        {
            arma::rowvec x_old = X.row(i);
            arma::mat x_new = mvrnormArma(1, x_old.t(), Sigma_prop);
            x_new = arma::conv_to<arma::rowvec>::from(x_new);
            arma::mat X_new = X;
            X_new.row(i) = x_new;
            //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& variatior a<1>
            X_new = standarize_matrix(X_new);
            x_new = X_new.row(i);
            //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            arma::vec curr_llkhd(1);
            curr_llkhd = log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                            X_new,
                                                            all_subset_mat,
                                                            param);
            arma::vec log_ratio(1);
            log_ratio = curr_llkhd;
            arma::mat x_new_in_mat = arma::conv_to<arma::mat>::from(x_new);
            arma::mat x_old_in_mat = arma::conv_to<arma::mat>::from(x_old);
            //log_ratio += dmvnrm_arma_fast(x_new_in_mat, x_old, Sigma_prop, true);   // proposal density // BUG
            log_ratio += dmvnrm_arma_fast(x_new_in_mat, mu_0, Sigma_0, true); // prior density of x_new in numerator  // CHEATING !!!
            log_ratio -= old_llkhd;
            //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& variation b<2>
            log_ratio -= dmvnrm_arma_fast(x_old_in_mat, mu_0, Sigma_0, true); // prior density of x_old in denominator // CHEATING !!!
            log_ratio -= dmvnrm_arma_fast(x_old_in_mat, x_old, Sigma_prop, true);
            //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            bool valid = all(log(arma::vec(1, arma::fill::randu)) < log_ratio);
            if (valid)
            {
                X = X_new;
                old_llkhd = curr_llkhd;
            }
            else
            {
                rej += 1;
            }
        }
        //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& variation c<3>
        Feat_sweep_lst[j] = standarize_matrix(X);
        //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    }
    int warm_start = ceil(m_t / 2);
    Rcpp::List Feat_sweep_lst_warm = Feat_sweep_lst[Rcpp::Range(warm_start, m_t - 1)];
    arma::vec likelihood_estep(1, arma::fill::zeros);
    for (int t = 0; t < Feat_sweep_lst.size(); t++)
    {
        likelihood_estep += log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                               Feat_sweep_lst[t],
                                                               all_subset_mat,
                                                               param);
    }
    likelihood_estep = likelihood_estep / Feat_sweep_lst.size();
    Rcpp::List res = Rcpp::List::create();
    res["Feat.lst"] = Feat_sweep_lst;
    res["likelihood.estep"] = likelihood_estep;
    res["last.feat.mat"] = X;
    res["rej"] = rej;
    return (res);
}

// ****************************************************************************************************************
// **************************************************************************************************************** Obsolete !!!
// Function estep_2 in mcem
// [[Rcpp::export]]
Rcpp::List mcem_estep_cpp_a1b2c2(const int em_round, const int min_sweeps,
                                 arma::mat Sigma_prop,
                                 Rcpp::List res_from_last_step,
                                 const int tot_author, const arma::imat all_subset_mat, const arma::irowvec data_s,
                                 bool if_message_reported)
{
    int rej = 0;
    Rcpp::List beta_lst = res_from_last_step["beta.lst"];
    int len_beta_lst = beta_lst.size();
    double beta_last = beta_lst[len_beta_lst - 1];
    arma::vec param(2, arma::fill::ones);
    param[1] = beta_last;
    arma::mat X = res_from_last_step["last.feat.mat"];
    arma::mat X_pivot = X;
    int num_sweeps = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List llkhd_mstep = res_from_last_step["likelihood.mstep.lst"];
    arma::vec old_llkhd(1);
    old_llkhd[0] = llkhd_mstep[em_round - 2];
    arma::rowvec mu_0 = arma::rowvec(2, arma::fill::zeros);
    arma::mat Sigma_0 = arma::mat(2, 2, arma::fill::eye);
    int m_t = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List Feat_sweep_lst(m_t);
    arma::mat X_anchor = res_from_last_step["anchor.feat.mat"];
    for (arma::uword j = 0; j < m_t; j++)
    {
        for (int i = 0; i < tot_author; i++)
        {
            arma::rowvec x_old = X.row(i);
            arma::mat x_new = mvrnormArma(1, x_old.t(), Sigma_prop);
            x_new = arma::conv_to<arma::rowvec>::from(x_new);
            arma::mat X_new = X;
            X_new.row(i) = x_new;
            if (false)
            { //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Tuning block
                // do nothing
            }
            else
            {
                X_new = standarize_matrix(X_new);
                x_new = X_new.row(i);
            } //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            arma::vec curr_llkhd(1);
            //Rcpp::print(Rcpp::wrap(data_s));
            curr_llkhd = log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                            X_new,
                                                            all_subset_mat,
                                                            param);
            arma::vec log_ratio(1);
            log_ratio = curr_llkhd;
            arma::mat x_new_in_mat = arma::conv_to<arma::mat>::from(x_new);
            arma::mat x_old_in_mat = arma::conv_to<arma::mat>::from(x_old);
            log_ratio += dmvnrm_arma_fast(x_new_in_mat, x_old, Sigma_prop, true); // proposal density
            log_ratio -= old_llkhd;
            if (true)
            {                                                                     //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Tuning block
                log_ratio -= dmvnrm_arma_fast(x_old_in_mat, mu_0, Sigma_0, true); // CHEATING !!! prior density
            }
            else
            {
                log_ratio -= dmvnrm_arma_fast(x_old_in_mat, x_old, Sigma_prop, true);
            } //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            bool valid = all(log(arma::vec(1, arma::fill::randu)) < log_ratio);
            if (valid)
            {
                X = X_new;
                old_llkhd = curr_llkhd;
            }
            else
            {
                rej += 1;
            }
        }
        if (true)
        {                          //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Tuning block
            Feat_sweep_lst[j] = X; // c<3>
        }
        else
        {
            arma::mat X_prcr = find_procrustean_transform(X_anchor, X);
            Feat_sweep_lst[j] = X_prcr;
        } //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    }
    int warm_start = ceil(m_t / 2);
    Rcpp::List Feat_sweep_lst_warm = Feat_sweep_lst[Rcpp::Range(warm_start, m_t - 1)];
    arma::vec likelihood_estep(1, arma::fill::zeros);
    for (int t = 0; t < Feat_sweep_lst.size(); t++)
    {
        likelihood_estep += log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                               Feat_sweep_lst[t],
                                                               all_subset_mat,
                                                               param);
    }
    likelihood_estep = likelihood_estep / Feat_sweep_lst.size();
    Rcpp::List res = Rcpp::List::create();
    res["Feat.lst"] = Feat_sweep_lst;
    res["likelihood.estep"] = likelihood_estep;
    res["last.feat.mat"] = X;
    res["rej"] = rej;
    return (res);
}

// Function estep_2 in mcem
// [[Rcpp::export]]
Rcpp::List mcem_estep_cpp_a1b2c4(const int em_round, const int min_sweeps,
                                 arma::mat Sigma_prop,
                                 Rcpp::List res_from_last_step,
                                 const int tot_author, const arma::imat all_subset_mat, const arma::irowvec data_s,
                                 bool if_message_reported)
{
    int rej = 0;
    Rcpp::List beta_lst = res_from_last_step["beta.lst"];
    int len_beta_lst = beta_lst.size();
    double beta_last = beta_lst[len_beta_lst - 1];
    arma::vec param(2, arma::fill::ones);
    param[1] = beta_last;
    arma::mat X = res_from_last_step["last.feat.mat"];
    arma::mat X_pivot = X;
    int num_sweeps = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List llkhd_mstep = res_from_last_step["likelihood.mstep.lst"];
    arma::vec old_llkhd(1);
    old_llkhd[0] = llkhd_mstep[em_round - 2];
    int m_t = max(Rcpp::IntegerVector::create(min_sweeps, pow(em_round, 2)));
    Rcpp::List Feat_sweep_lst(m_t);
    arma::mat X_anchor = res_from_last_step["anchor.feat.mat"];
    for (arma::uword j = 0; j < m_t; j++)
    {
        for (int i = 0; i < tot_author; i++)
        {
            arma::rowvec x_old = X.row(i);
            arma::mat x_new = mvrnormArma(1, x_old.t(), Sigma_prop);
            x_new = arma::conv_to<arma::rowvec>::from(x_new);
            arma::mat X_new = X;
            X_new.row(i) = x_new;
            if (false)
            { //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Tuning block
                // do nothing
            }
            else
            {
                X_new = standarize_matrix(X_new);
                x_new = X_new.row(i);
            } //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            arma::vec curr_llkhd(1);
            //Rcpp::print(Rcpp::wrap(data_s));
            curr_llkhd = log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                            X_new,
                                                            all_subset_mat,
                                                            param);
            arma::vec log_ratio(1);
            log_ratio = curr_llkhd;
            arma::mat x_new_in_mat = arma::conv_to<arma::mat>::from(x_new);
            arma::mat x_old_in_mat = arma::conv_to<arma::mat>::from(x_old);
            log_ratio += dmvnrm_arma_fast(x_new_in_mat, x_old, Sigma_prop, true);
            log_ratio -= old_llkhd;
            if (true)
            { //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Tuning block
                //log_ratio -= dmvnrm_arma_fast(x_old_in_mat, arma::rowvec(2, arma::fill::zeros), Sigma_prop, true);
                log_ratio -= dmvnrm_arma_fast(x_old_in_mat, arma::rowvec(2, arma::fill::zeros), arma::mat(2, 2, arma::fill::eye), true);
                //Rcpp::print(Rcpp::wrap(dmvnrm_arma_fast(x_old_in_mat, arma::rowvec(2, arma::fill::zeros), arma::mat(2,2,arma::fill::eye), true)));
                //Rcpp::print(Rcpp::wrap(dmvnrm_arma_fast(x_old_in_mat, x_old, Sigma_prop, true)));
            }
            else
            {
                log_ratio -= dmvnrm_arma_fast(x_old_in_mat, x_old, Sigma_prop, true);
            } //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            bool valid = all(log(arma::vec(1, arma::fill::randu)) < log_ratio);
            if (valid)
            {
                X = X_new;
                old_llkhd = curr_llkhd;
            }
            else
            {
                rej += 1;
            }
        }
        if (true)
        { //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Tuning block
            arma::mat X_prcr = find_procrustean_transform(X_anchor, X);
            Feat_sweep_lst[j] = X_prcr;
        }
        else
        {
            arma::mat X_prcr = find_procrustean_transform(X_anchor, X);
            Feat_sweep_lst[j] = X_prcr;
        } //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    }
    int warm_start = ceil(m_t / 2);
    Rcpp::List Feat_sweep_lst_warm = Feat_sweep_lst[Rcpp::Range(warm_start, m_t - 1)];
    arma::vec likelihood_estep(1, arma::fill::zeros);
    for (int t = 0; t < Feat_sweep_lst.size(); t++)
    {
        likelihood_estep += log_likelihood_shape_cvxhull_gamma(data_s.t(),
                                                               Feat_sweep_lst[t],
                                                               all_subset_mat,
                                                               param);
    }
    likelihood_estep = likelihood_estep / Feat_sweep_lst.size();
    Rcpp::List res = Rcpp::List::create();
    res["Feat.lst"] = Feat_sweep_lst;
    res["likelihood.estep"] = likelihood_estep;
    res["last.feat.mat"] = X;
    res["rej"] = rej;
    return (res);
}
