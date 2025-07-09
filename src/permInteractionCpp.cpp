#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat permInteractionCpp(
  const arma::mat& XY,           // Note: Use const ref for input, unless you intend to modify
  const arma::ivec& g_XY,        // ditto, and safe for Armadillo subsetting
  int n1,
  int n2,
  int B,
  int g1,
  int g2
) {

  int n = XY.n_rows;
  int p = XY.n_cols;
  arma::mat T_perm(B, p, fill::zeros);

  arma::ivec idx = regspace<arma::ivec>(0, n - 1);

  for (int b = 0; b < B; ++b) {
    arma::ivec perm_idx = shuffle(idx);
    arma::uvec p1_u = arma::conv_to<arma::uvec>::from(perm_idx.subvec(0, n1 - 1));
    arma::uvec p2_u = arma::conv_to<arma::uvec>::from(perm_idx.subvec(n1, n1 + n2 - 1));

    // Subset g_XY for each group (the genotype vector for samples in group)
    arma::ivec g1_p1 = g_XY.elem(p1_u);
    arma::ivec g1_p2 = g_XY.elem(p2_u);

    arma::uvec mask_g1_p1 = find(g1_p1 == g1);
    arma::uvec mask_g2_p1 = find(g1_p1 == g2);
    arma::uvec mask_g1_p2 = find(g1_p2 == g1);
    arma::uvec mask_g2_p2 = find(g1_p2 == g2);

    if (mask_g1_p1.is_empty() || mask_g2_p1.is_empty() ||
        mask_g1_p2.is_empty() || mask_g2_p2.is_empty()) {
      T_perm.row(b).fill(datum::nan);
      continue;
    }

    arma::mat XY_p1 = XY.rows(p1_u);
    arma::mat XY_p2 = XY.rows(p2_u);

    arma::rowvec m_g1_p1 = mean(XY_p1.rows(mask_g1_p1), 0);
    arma::rowvec m_g2_p1 = mean(XY_p1.rows(mask_g2_p1), 0);
    arma::rowvec m_g1_p2 = mean(XY_p2.rows(mask_g1_p2), 0);
    arma::rowvec m_g2_p2 = mean(XY_p2.rows(mask_g2_p2), 0);

    arma::rowvec dX = m_g2_p1 - m_g1_p1;
    arma::rowvec dY = m_g2_p2 - m_g1_p2;

    T_perm.row(b) = dY - dX;
  }
  return T_perm;
}
