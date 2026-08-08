#pragma once

#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <cstdint>
#include <string>
#include <vector>
#include <limits>

enum EnumStateTrackStrategy
{
  ACTIVE_SET_CHANGE = 0,
  ALL_ITERATION = 1,
  EVERY_K_ITERATION = 2,
  NO_STATE_TRACKING = 3
};

enum EnumStagewisePhase
{
  STAGEWISE_NOT_STARTED = 0,
  STAGEWISE_NONPEN = 1,
  STAGEWISE_SATURATED = 2,
  STAGEWISE_ITERATION = 3,
  STAGEWISE_FINISHED = 4
};
enum EnumStagewiseTerminationReason
{
  STAGEWISE_NOT_INITIALIZED = 0,
  STAGEWISE_RUNNING = 1,

  STAGEWISE_BETA_STEP_ZERO = 2,
  STAGEWISE_BETA_STALLED = 3,
  STAGEWISE_OBJECTIVE_STALLED = 4,
  STAGEWISE_PSEUDO_R2_CUTOFF_REACHED = 5,

  STAGEWISE_EPSILON_MIN_REACHED = 6,
  STAGEWISE_ITERATION_LIMIT_REACHED = 7
};

struct IEgmifsCriterion;
struct EgmifsInput;
struct EgmifsControl;

struct IEgmifsLinkFunc
{
  virtual ~IEgmifsLinkFunc() = default;

  virtual std::string name() const = 0;

  virtual arma::uword parameter_count() const noexcept = 0;

  virtual void prepare(
      const EgmifsInput& input,
      const EgmifsControl& control
  ) const = 0;

  virtual void inverse(
      const arma::vec& eta,
      const arma::vec& link_parameters,
      arma::vec& mu
  ) const = 0;

  virtual void grad(
      const arma::vec& eta,
      const arma::vec& link_parameters,
      arma::vec& d_mu_d_eta,
      arma::mat& d_mu_d_link_parameters
  ) const = 0;

  virtual arma::vec initial_parameters() const
  {
    return arma::vec(parameter_count(), arma::fill::zeros);
  }

  virtual arma::vec parameter_lower_bounds() const
  {
    arma::vec lower(parameter_count());
    lower.fill(-arma::datum::inf);
    return lower;
  }

  virtual arma::vec parameter_upper_bounds() const
  {
    arma::vec upper(parameter_count());
    upper.fill(arma::datum::inf);
    return upper;
  }
};
struct IEgmifsFamily
{
  virtual ~IEgmifsFamily() = default;

  virtual std::string name() const = 0;

  virtual arma::uword parameter_count() const noexcept = 0;

  virtual void prepare(
      const EgmifsInput& input,
      const EgmifsControl& control
  ) const = 0;

  virtual void negloglik(
      const arma::vec& y,
      const arma::vec& mu,
      const arma::vec& family_parameters,
      double& negloglik
  ) const = 0;

  virtual void grad(
      const arma::vec& y,
      const arma::vec& mu,
      const arma::vec& family_parameters,
      arma::vec& d_negloglik_d_mu,
      arma::vec& d_negloglik_d_family_parameters
  ) const = 0;

  virtual arma::vec initial_parameters() const
  {
    return arma::vec(parameter_count(), arma::fill::zeros);
  }

  virtual arma::vec parameter_lower_bounds() const
  {
    arma::vec lower(parameter_count());
    lower.fill(-arma::datum::inf);
    return lower;
  }

  virtual arma::vec parameter_upper_bounds() const
  {
    arma::vec upper(parameter_count());
    upper.fill(arma::datum::inf);
    return upper;
  }
};


/*
 * Combined family-link interface.
 *
 * A custom implementation may evaluate d(negative log-likelihood)/d(eta)
 * directly using a numerically stable fused formula. When no combined
 * implementation is supplied, the package constructs an internal adapter
 * that delegates to IEgmifsFamily and IEgmifsLinkFunc and applies
 * the ordinary chain rule.
 */
struct IEgmifsFamilyLink
{
  virtual ~IEgmifsFamilyLink() = default;

  virtual std::string family_name() const = 0;
  virtual std::string link_name() const = 0;

  virtual arma::uword family_parameter_count() const noexcept = 0;
  virtual arma::uword link_parameter_count() const noexcept = 0;

  virtual arma::vec family_initial_parameters() const
  {
    return arma::vec(
      family_parameter_count(),
      arma::fill::zeros
    );
  }

  virtual arma::vec family_parameter_lower_bounds() const
  {
    arma::vec lower(family_parameter_count());
    lower.fill(-arma::datum::inf);
    return lower;
  }

  virtual arma::vec family_parameter_upper_bounds() const
  {
    arma::vec upper(family_parameter_count());
    upper.fill(arma::datum::inf);
    return upper;
  }

  virtual arma::vec link_initial_parameters() const
  {
    return arma::vec(
      link_parameter_count(),
      arma::fill::zeros
    );
  }

  virtual arma::vec link_parameter_lower_bounds() const
  {
    arma::vec lower(link_parameter_count());
    lower.fill(-arma::datum::inf);
    return lower;
  }

  virtual arma::vec link_parameter_upper_bounds() const
  {
    arma::vec upper(link_parameter_count());
    upper.fill(arma::datum::inf);
    return upper;
  }

  virtual void prepare(
      const EgmifsInput& input,
      const EgmifsControl& control
  ) const = 0;

  virtual void inverse(
      const arma::vec& eta,
      const arma::vec& link_parameters,
      arma::vec& mu
  ) const = 0;

  virtual void negloglik(
      const arma::vec& y,
      const arma::vec& mu,
      const arma::vec& family_parameters,
      double& negloglik
  ) const = 0;

  /*
   * Fill all derivatives needed by the fitting code.
   *
   * d_negloglik_d_eta is the derivative used for beta and theta. A custom
   * combined implementation may calculate it directly rather than as
   * d_negloglik_d_mu % d_mu_d_eta.
   */
  virtual void grad(
      const arma::vec& y,
      const arma::vec& eta,
      const arma::vec& mu,
      const arma::vec& family_parameters,
      const arma::vec& link_parameters,
      arma::vec& d_negloglik_d_mu,
      arma::vec& d_mu_d_eta,
      arma::mat& d_mu_d_link_parameters,
      arma::vec& d_negloglik_d_eta,
      arma::vec& d_negloglik_d_family_parameters,
      arma::vec& d_negloglik_d_link_parameters
  ) const = 0;
};



/*
 * Public read-only input struct.
 *
 * This exposes model data to user-defined C++ criteria.
 * The matrix/vector fields are const references, so plugins can read but should
 * not mutate package-owned data.
 */
struct EgmifsInput {
  const arma::mat& X;
  const arma::vec& y;
  const arma::mat& w;
  const arma::vec& offset;
  const arma::vec& weight_vec;
  bool has_prior;
  double enet_alpha;

  /*
   * family_link is always non-null internally.
   *
   * If the user supplies a combined implementation, family and link_func may
   * be null. Otherwise family/link_func are retained and family_link points to
   * an internal delegating adapter.
   */
  const IEgmifsFamily* family;
  const IEgmifsLinkFunc* link_func;
  const IEgmifsFamilyLink* family_link;
  std::vector<const IEgmifsCriterion*> criteria;
};


/*
 * Public read-only control view.
 *
 * This exposes model data to user-defined C++ criteria.
 * The matrix/vector fields are const references, so plugins can read but should
 * not mutate package-owned data.
 */
struct EgmifsNloptControl
{
  nlopt_algorithm algorithm;
  double xtol_rel;
  double ftol_rel;
  int maxeval;
};

struct EgmifsControl {
  uint64_t null_iteration_max;
  uint64_t stagewise_iteration_max;

  double null_family_parameter_abs_tol;
  double stagewise_objective_rel_tol;
  double stagewise_beta_step_norm_tol;

  double epsilon_max;
  double epsilon_start;
  double epsilon_min;
  double loglik_reltol_cutoff;

  double enet_abs_tol;
  double enet_rel_tol;
  uint32_t enet_max_iter;

  EnumStateTrackStrategy state_track_strategy;
  uint64_t state_track_freq;
  bool verbose;
  bool include_data;

  const arma::vec& theta_initial;
  const arma::vec& theta_lower_bounds;
  const arma::vec& theta_upper_bounds;

  EgmifsNloptControl nonpen_nlopt;
  EgmifsNloptControl family_nlopt;
  EgmifsNloptControl link_nlopt;
};

struct EgmifsParameters
{
  arma::vec beta;
  arma::vec theta;

  // Future module-specific parameters passed from loglik/other parts
  arma::vec family_parameters;
  arma::vec link_parameters;
};

struct EgmifsGradients
{
  arma::vec d_negloglik_d_mu;
  arma::vec d_mu_d_eta;
  arma::vec d_negloglik_d_eta;

  const arma::mat& d_eta_d_beta;
  const arma::mat& d_eta_d_theta;

  // Future gradients wrt other parameters such as dispersion
  arma::mat d_mu_d_link_parameters;
  arma::vec d_negloglik_d_family_parameters;
  arma::vec d_negloglik_d_link_parameters;

  // main gradient for FS
  arma::vec d_negloglik_d_beta;
};

struct EgmifsPredictors
{
  EgmifsParameters param;

  arma::vec xbeta;
  arma::vec wtheta;
  arma::vec eta;
  arma::vec mu;

  arma::uvec active_set;
};

struct EgmifsStagewise
{
  EnumStagewisePhase phase = EnumStagewisePhase::STAGEWISE_NOT_STARTED;
  EnumStagewiseTerminationReason termination_reason =
    EnumStagewiseTerminationReason::STAGEWISE_NOT_INITIALIZED;

  /*
   * Immutable beta at the beginning of the current
   * stagewise iteration.
   */
  arma::vec beta_start;

  /*
   * Candidate beta for the current epsilon:
   *
   *   beta_trial = beta_start + delta_beta
   */
  arma::vec beta_trial;

  /*
   * Elastic-net step for the current epsilon.
   */
  arma::vec delta_beta;

  /*
   * Generic fitting workspaces retained across optimizer calls.
   * These belong to the stagewise state, while optimizer objects and
   * callback-specific bookkeeping remain implementation details.
   */
  arma::vec saturated_family_parameters;
  double saturated_negloglik = arma::datum::nan;

  arma::vec theta_before_optimize;
  arma::vec family_parameters_before_optimize;
  arma::vec link_parameters_before_optimize;
  arma::vec saturated_family_parameters_before_optimize;

  /*
   * Numerical explanation produced by stagewise.
   *
   * Examples:
   *   "Relative negative log-likelihood difference = ..."
   *   "Epsilon = ..., epsilon_min = ..."
   */
  std::string termination_detail = "Not initialized";

  double epsilon = 0.0;

  double negloglik_trial = arma::datum::nan;

  uint32_t halving_count = 0;
};

struct EgmifsState
{
  EgmifsPredictors param;

  double negloglik = arma::datum::nan;

  Rcpp::NumericVector criteria;

  uint64_t iteration = 0;

  /*
   * This value is populated when a State snapshot is saved to Path.
   * The live mutable State does not own the null/saturated baselines.
   */
  double pseudo_r2 = arma::datum::nan;

  /*
   * Wall-clock time spent completing this stagewise iteration, in seconds.
   * The null-model snapshot at iteration 0 keeps this value at zero because
   * null-model timing is stored once at Path level.
   */
  double elapsed_time = 0.0;
};

struct EgmifsBestCriterion
{
  std::string name;

  double value =
    std::numeric_limits<double>::infinity();

  /*
   * Zero-based index into the R-backed saved-state pool in EgmifsPath::saved_states.
   * A negative value means that the best state is still pending and has not
   * yet been materialized.
   */
  R_xlen_t state_index = -1;
};

struct EgmifsSavedStates
{
  /*
   * Capacity-backed R storage for every distinct materialized iteration.
   * Ordinary path tracking and best criteria refer to these entries by
   * zero-based index.
   */
  Rcpp::CharacterVector state_names;
  Rcpp::NumericVector iterations;
  Rcpp::NumericVector negloglik;
  Rcpp::NumericVector pseudo_r2;
  Rcpp::NumericVector elapsed_time;

  Rcpp::List criteria;
  Rcpp::List beta;
  Rcpp::List theta;
  Rcpp::List family_parameters;
  Rcpp::List link_parameters;
  Rcpp::List xbeta;
  Rcpp::List wtheta;
  Rcpp::List eta;
  Rcpp::List mu;
  Rcpp::List active_set;

  /* Number of initialized entries in the capacity-backed storage. */
  R_xlen_t count = 0;
};

struct EgmifsPath
{
  /*
   * Fit-wide null-model result. Beta is omitted because the null model
   * always has beta = 0 and storing a p-length zero vector is unnecessary.
   */
  double null_negloglik = arma::datum::nan;
  arma::vec null_theta;
  arma::vec null_family_parameters;
  arma::vec null_link_parameters;

  /*
   * Fit-wide saturated-model result. Family parameters are kept generic:
   * an empty vector for parameter-free families, one value for NB2, and
   * any other size required by another family implementation.
   */
  double saturated_negloglik = arma::datum::nan;
  arma::vec saturated_family_parameters;

  /* Wall-clock fitting times, in seconds. */
  double null_time = 0.0;
  double saturated_time = 0.0;
  double total_time = 0.0;

  std::vector<EgmifsBestCriterion>
    best_criteria;

  /*
   * Every distinct saved iteration is materialized at most once here.
   * Ordinary path tracking and all best criteria refer to it by index.
   */
  EgmifsSavedStates saved_states;

  /* Zero-based indices selecting the ordinary tracked path states. */
  Rcpp::IntegerVector state_indices;

  /* Number of initialized entries in state_indices. */
  R_xlen_t state_count = 0;

  arma::uvec last_saved_active_set;
  bool active_set_changed = false;
  std::string message;
};

struct EgmifsRuntime {
  const EgmifsState& state;
  const EgmifsPath& path;
  const EgmifsGradients& gradient;
  const EgmifsStagewise& stagewise;
};

struct EgmifsContext {
  const EgmifsInput& input;
  const EgmifsControl& control;
  const EgmifsRuntime& runtime;
};

struct IEgmifsCriterion
{
  virtual ~IEgmifsCriterion() = default;

  virtual std::string name() const = 0;

  virtual void prepare(
      const EgmifsInput& input,
      const EgmifsControl& control
  ) const = 0;

  virtual double evaluate(
      const EgmifsInput& input,
      const EgmifsControl& control,
      const EgmifsState& state
  ) const = 0;
};
