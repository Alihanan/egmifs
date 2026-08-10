#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]

#include <cstdint>

#include "context.h"


// [[Rcpp::export]]
Rcpp::List egmifs_cpp(
    const arma::mat& X,
    const arma::vec& y,
    const arma::mat& w,
    const arma::vec& offset,

    const arma::vec& weight_vec,
    double enet_alpha,

    double epsilon_start,
    double epsilon_max,
    double epsilon_min,

    uint32_t null_iteration_max,
    uint32_t stagewise_iteration_max,
    double null_family_parameter_abs_tol,
    double stagewise_objective_rel_tol,
    double stagewise_beta_step_norm_tol,

    SEXP family,
    SEXP link_func,
    Rcpp::Nullable<Rcpp::List> criteria,

    double loglik_reltol_cutoff,
    double enet_abs_tol,
    double enet_rel_tol,
    uint32_t enet_max_iter,

    bool verbose,
    int state_track_strategy,
    uint64_t state_track_freq,
    bool include_data,

    const arma::vec& theta_initial,
    const arma::vec& theta_lower_bounds,
    const arma::vec& theta_upper_bounds,

    int nonpen_nlopt_algorithm,
    double nonpen_nlopt_xtol_rel,
    double nonpen_nlopt_ftol_rel,
    int nonpen_nlopt_maxeval,

    int family_nlopt_algorithm,
    double family_nlopt_xtol_rel,
    double family_nlopt_ftol_rel,
    int family_nlopt_maxeval,

    int link_nlopt_algorithm,
    double link_nlopt_xtol_rel,
    double link_nlopt_ftol_rel,
    int link_nlopt_maxeval,

    SEXP family_link = R_NilValue
) {
  EgmifsContextInternal ctx(
    X,
    y,
    w,
    offset,

    weight_vec,
    enet_alpha,
    epsilon_start,
    epsilon_max,
    epsilon_min,
    null_iteration_max,
    stagewise_iteration_max,
    null_family_parameter_abs_tol,
    stagewise_objective_rel_tol,
    stagewise_beta_step_norm_tol,

    family,
    link_func,
    criteria,

    loglik_reltol_cutoff,
    enet_abs_tol,
    enet_rel_tol,
    enet_max_iter,

    verbose,
    state_track_strategy,
    state_track_freq,
    include_data,

    theta_initial,
    theta_lower_bounds,
    theta_upper_bounds,

    nonpen_nlopt_algorithm,
    nonpen_nlopt_xtol_rel,
    nonpen_nlopt_ftol_rel,
    nonpen_nlopt_maxeval,

    family_nlopt_algorithm,
    family_nlopt_xtol_rel,
    family_nlopt_ftol_rel,
    family_nlopt_maxeval,

    link_nlopt_algorithm,
    link_nlopt_xtol_rel,
    link_nlopt_ftol_rel,
    link_nlopt_maxeval,

    family_link
  );

  ctx.stagewise.fit();

  Rcpp::List output = ctx.to_list();

  output.attr("class") =
    Rcpp::CharacterVector::create(
      "egmifs",
      "list"
    );

  return output;
}
