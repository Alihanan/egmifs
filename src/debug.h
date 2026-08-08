#pragma once

#define EGMIFS_VERBOSE(VERBOSE, MSG)                          \
do {                                                               \
  if (VERBOSE) {                                                   \
    Rcpp::Rcout << "[egmifs] " << MSG << "\n";                \
  }                                                                \
} while (false)
