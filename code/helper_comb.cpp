// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//#include <iostream>

// Function 0
//      A recursive function helper used in find_subset(), update output matrix by appending valid k-sets row by row.
//      adopted from https://www.sanfoundry.com/cpp-program-generate-subsets-k-length/
// [[Rcpp::export]]
void find_subset_recursive(arma::imat &output, arma::irowvec arr, int n, int k, int start_pos, int currt_pos, arma::uvec if_used)
{
    if (currt_pos > k)
        return;
    else if (currt_pos == k)
    {
        arma::irowvec k_set(k);
        int count = 0;
        for (int i = 0; i < n; i++)
        {
            if (if_used[i] == true)
            {
                k_set[count] = arr[i];
                count += 1;
            }
        }
        output = join_vert(output, k_set);
        return;
    }
    if (start_pos == n)
    {
        return;
    }
    if_used[start_pos] = true;
    find_subset_recursive(output, arr, n, k, start_pos + 1, currt_pos + 1, if_used);
    if_used[start_pos] = false;
    find_subset_recursive(output, arr, n, k, start_pos + 1, currt_pos, if_used);
}

// Function 1: Find all k-sets from {1,2,...,n}.
// [[Rcpp::export]]
arma::imat find_subset(int n, int k)
{
    arma::uvec if_used(n);
    arma::irowvec arr(n);
    for (int i = 0; i < n; i++)
    {
        arr[i] = i + 1;
        if_used[i] = false;
    }
    //std::cout << "\nFind all " << k << "-sets out of 1 :" << n << " using Rcpp \n ";
    arma::imat output(1, k, arma::fill::zeros);
    find_subset_recursive(output, arr, n, k, 0, 0, if_used);
    output.shed_row(0);
    return output;
}
