library(testthat);
library(Rssa);

compute.reconstructions <- function(x, Ls, groups,
                                    kind = c("1d-ssa", "toeplitz-ssa"),
                                    svd.method = c("eigen", "propack", "nutrlan", "svd"),
                                    neig = max(unlist(groups)) + 1) {
  kind <- match.arg(kind);
  svd.method <- match.arg(svd.method);
  if (is.null(neig))
    neig <- eval(formals()$neig);

  out <- lapply(Ls, function(L) {
        ss <- if (identical(svd.method, "eigen") || identical(svd.method, "svd")) {
              ssa(x, L, kind = kind, svd.method = svd.method)
            } else {
              suppressWarnings(ssa(x, L, kind = kind, svd.method = svd.method, neig = neig));
            }
        reconstruct(ss, groups = groups);
      });

  names(out) <- paste("L", Ls, sep = "");

  attributes(out)[c("kind")] <- list(kind);
  out;
}

compute.forecasts <- function(x, Ls, groups, len,
                              kind = c("1d-ssa", "toeplitz-ssa"),
                              forecast.method = c("recurent", "vector"),
                              base = c("reconstructed", "original"),
                              svd.method = c("eigen", "propack", "nutrlan", "svd"),
                              neig = max(unlist(groups)) + 1) {
  kind <- match.arg(kind);
  svd.method <- match.arg(svd.method);
  forecast.method <- match.arg(forecast.method);
  base <- match.arg(base);
  if (is.null(neig))
    neig <- eval(formals()$neig);

  out <- lapply(Ls, function(L) {
        ss <- if (identical(svd.method, "eigen") || identical(svd.method, "svd")) {
              ssa(x, L, kind = kind, svd.method = svd.method)
            } else {
              suppressWarnings(ssa(x, L, kind = kind, svd.method = svd.method, neig = neig));
            }

        rec <- if (identical(forecast.method, "recurent")) {
              rforecast(ss, groups = groups, len = len, base = base, only.new = FALSE);
            } else if (identical(forecast.method, "vector")) {
              vforecast(ss, groups = groups, len = len, only.new = FALSE);
            }

        rec;
      });

  names(out) <- paste("L", Ls, sep = "");
  attributes(out)[c("kind", "forecast.method", "base")] <- list(kind, forecast.method, base);
  out;
}

make.test.data <- function(what = c("reconstruct", "rforecast", "vforecast"),
                           series,
                           name = deparse(substitute(series)),
                           Ls,
                           Ls.forecast = Ls,
                           groups,
                           groups.forecast = groups,
                           len = 100,
                           kind = c("1d-ssa", "toeplitz-ssa"),
                           svd.method = c("eigen", "propack", "svd", "nutrlan"),
                           neig = NULL,
                           tolerance = 1e-7,
                           svd.methods = c("eigen", "propack", "svd", "nutrlan")) {
  what <- sapply(what, match.arg, choices = eval(formals()$what));
  kind <- match.arg(kind);
  svd.method <- match.arg(svd.method);
  svd.methods <- sapply(svd.methods, match.arg, choices = eval(formals()$svd.methods));

  out <- list(series = series);

  if ("reconstruct" %in% what) {
    out$reconstruction <- compute.reconstructions(series, Ls, groups = groups,
                                                  kind,
                                                  svd.method = svd.method, neig = neig);
  }

  if ("rforecast" %in% what) {
    out$rforecast.orig <- compute.forecasts(series, Ls.forecast, groups = groups.forecast, len = len,
                                            kind = kind, forecast.method = "recurent", base = "original",
                                            svd.method = svd.method, neig = neig);
    out$rforecast.rec <- compute.forecasts(series, Ls.forecast, groups = groups.forecast, len = len,
                                           kind = kind, forecast.method = "recurent", base = "reconstructed",
                                           svd.method = svd.method, neig = neig);
  }

  if ("vforecast" %in% what) {
    out$vforecast <- compute.forecasts(series, Ls.forecast, groups = groups.forecast, len = len,
                                       kind = kind, forecast.method = "vector",
                                       svd.method = svd.method, neig = neig);
  }

  attr(out, "name") <- name;
  attr(out, "what") <- what;
  attr(out, "pars") <- list(kind = kind,
                            Ls = Ls,
                            groups = groups,
                            Ls.forecast = Ls.forecast,
                            groups.forecast = groups.forecast,
                            len = len,
                            neig = neig);
  attr(out, "tolerance") <- tolerance;
  attr(out, "svd.methods") <- svd.methods;

  out;
}

test.test.data <- function(what,
                           test.data,
                           name = attr(test.data, "name"),
                           svd.methods,
                           neig,
                           tolerance,
                           ...) {
  if (missing(tolerance))
    tolerance <- attr(test.data, "tolerance");

  pars <- attr(test.data, "pars");

  if (missing(neig))
    neig <- pars$neig;

  kind <- pars$kind;

  if (missing(what)) {
    what <- attr(test.data, "what");
  } else {
    what <- sapply(what, match.arg, choices = c("reconstruct", "rforecast", "vforecast"));
    what <- intersect(what, attr(test.data, "what"));
  }

  if (missing(svd.methods)) {
    svd.methods <- attr(test.data, "svd.methods");
  } else {
    svd.methods <- sapply(svd.methods, match.arg, choices = c("eigen", "propack", "svd", "nutrlan"));
  }

  series <- test.data$series;
  Ls <- pars$Ls;
  groups <- pars$groups;
  groups.forecast <- pars$groups.forecast;
  Ls.forecast <- pars$Ls.forecast;
  len <- pars$len;

  for (svd.method in svd.methods) {
    if ("reconstruct" %in% what) {
      reconstruction <- compute.reconstructions(series, Ls, groups = groups,
                                                kind,
                                                svd.method = svd.method, neig = neig);

      expect_equal(reconstruction, test.data$reconstruction,
                   label = sprintf("%s, %s: %s$reconstruction", name, kind, svd.method),
                   tolerance = tolerance, ...);
    }

    if ("rforecast" %in% what) {
      rforecast.orig <- compute.forecasts(series, Ls.forecast, groups = groups.forecast, len = len,
                                          kind = kind, forecast.method = "recurent", base = "original",
                                          svd.method = svd.method, neig = neig);
      rforecast.rec <- compute.forecasts(series, Ls.forecast, groups = groups.forecast, len = len,
                                         kind = kind, forecast.method = "recurent", base = "reconstructed",
                                         svd.method = svd.method, neig = neig);

      expect_equal(rforecast.orig, test.data$rforecast.orig,
                   label = sprintf("%s, %s: %s$rforecast.orig", name, kind, svd.method),
                   tolerance = tolerance, ...);
      expect_equal(rforecast.rec, test.data$rforecast.rec,
                   label = sprintf("%s, %s: %s$rforecast.rec", name, kind, svd.method),
                   tolerance = tolerance, ...);
    }

    if ("vforecast" %in% what) {
      vforecast <- compute.forecasts(series, Ls.forecast, groups = groups.forecast, len = len,
                                     kind = kind, forecast.method = "vector",
                                     svd.method = svd.method, neig = neig);

      expect_equal(vforecast, test.data$vforecast,
                   label = sprintf("%s, %s: %s$vforecast", name, kind, svd.method),
                   tolerance = tolerance, ...);
    }
  }
}

is_multisets_approx_equal <- function(mset1, mset2, tol = .Machine$double.eps^0.5) {
  if (length(mset1) != length(mset2))
    return(FALSE);

  for (el in mset1) {
    i <- which.min(abs(mset2 - el));
    if (abs(mset2[i] - el) > tol)
      return(FALSE);
    mset2 <- mset2[-i];
  }

  return(length(mset2) == 0);
}
