library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("Signal parameters' estimation");

# test_that("parestimate.esprit works correctly for polynomial trends", {
# for (d in 0:3) {
#   N <- 40;
#   h <- t(hankel((1:N) ^ d, d + 1));
#   U <- svd(h)$u
#   par <- parestimate.esprit(U);
#
#   mu <- par$moduli * exp(pi * 2i / par$periods);
#   # Sad but true. Multiple eigenroots are evaluated with bad precision
#   expect_true(is_multisets_approx_equal(mu, rep(1, d + 1), tol = 0.001),
#               label = sprintf("Estimated characteristical roots for (1:N)^%d are correct", d));
#   }
# });

test_that("parestimate.esprit works correctly for two sines", {
  r <- 4;
  N <- 40;
  T1 <- 6;
  T2 <- 11;
  v <- sin(2 * pi * (1:N) / T1) + sin(2 * pi * (1:N) / T2);
  h <- t(hankel(v, r));
  U <- qr.Q(qr(h));
  par <- parestimate.esprit(U);
  mu <- par$moduli * exp(pi * 2i / par$periods);
  expectred.mu <- exp(pi * 2i / c(T1, -T1, T2, -T2));

  expect_true(is_multisets_approx_equal(mu, expectred.mu),
              label = sprintf("Est. ch. roots for sum of sines with periods %3.1f and %3.1f are correct", T1, T2));
})

test_that("parestimate works correctly for two sines", {
  r <- 4
  N <- 40
  T1 <- 6
  T2 <- 11
  v <- sin(2 * pi * (1:N) / T1) + sin(2 * pi * (1:N) / T2)
  ss <- ssa(v)
  for (solve.method in c("ls", "tls")) {
    par <- parestimate(ss, groups = list(sines = 1:4),
                       solve.method = solve.method)
    mu <- par$moduli * exp(pi * 2i / par$periods)
    expectred.mu <- exp(pi * 2i / c(T1, -T1, T2, -T2))

    expect_true(is_multisets_approx_equal(mu, expectred.mu),
                label = sprintf("Est. ch. roots for sum of sines with periods %3.1f and %3.1f are correct", T1, T2))
  }
})

test_that("parestimate works correctly for two sines in shaped case", {
  r <- 4
  N <- 80
  L <- 15
  T1 <- 6
  T2 <- 11
  v <- sin(2 * pi * (1:N) / T1) + sin(2 * pi * (1:N) / T2)
  v[20:25] <- NA
  v[50:55] <- NA

  ss <- ssa(v, L = L)
  for (solve.method in c("ls", "tls")) {
    par <- parestimate(ss, groups = list(sines = 1:4),
                       solve.method = solve.method)
    mu <- par$moduli * exp(pi * 2i / par$periods)
    expectred.mu <- exp(pi * 2i / c(T1, -T1, T2, -T2))

    expect_true(is_multisets_approx_equal(mu, expectred.mu),
                label = sprintf("Est. ch. roots for sum of sines with periods %3.1f and %3.1f are correct", T1, T2))
  }
})

test_that("parestimate works correctly for two sines in circular case", {
  r <- 4
  N <- 40
  T1 <- 5
  T2 <- 8
  v <- sin(2 * pi * (1:N) / T1) + sin(2 * pi * (1:N) / T2)
  ss <- ssa(v, L = N, circular = TRUE)
  for (solve.method in c("ls", "tls")) {
    par <- parestimate(ss,
                       groups = list(sines = 1:4),
                       solve.method = solve.method)
    mu <- par$roots
    expectred.mu <- exp(pi * 2i / c(T1, -T1, T2, -T2))

    expect_true(is_multisets_approx_equal(mu, expectred.mu),
                label = sprintf("Est. ch. roots for sum of sines with periods %3.1f and %3.1f are correct", T1, T2))
  }
})

test_that("parestimate works correctly for two sines in 2d case", {
  r <- 4
  N1 <- 55
  N2 <- 40
  T1 <- 5
  T2 <- 8
  v1 <- sin(2 * pi * (1:N1) / T1)
  v2 <- sin(2 * pi * (1:N2) / T2)
  mx <- outer(v1, v2)
  ss <- ssa(mx, kind = "2d-ssa")
  for (solve.method in c("ls", "tls")) {
    for (pairing.method in c("memp", "diag")) {
      par <- parestimate(ss, groups = list(sines = 1:4),
                         solve.method = solve.method,
                         pairing.method = pairing.method)
      lm <- par[[1]]$roots
      mu <- par[[2]]$roots
      expectred.lm <- exp(pi * 2i / c(T1, -T1, T1, -T1))
      expectred.mu <- exp(pi * 2i / c(T2, -T2, T2, -T2))

      expect_true(is_multisets_approx_equal(mu, expectred.mu) && is_multisets_approx_equal(lm , expectred.lm),
                  label = sprintf("Est. ch. roots for sum of sines with periods %3.1f and %3.1f are correct (method = %s-%s)",
                                  T1, T2, solve.method, pairing.method))
    }
  }
})

test_that("parestimate works correctly for two sines in circular 2d case", {
  r <- 4
  N1 <- 55
  N2 <- 40
  T1 <- 5
  T2 <- 8
  expectred.lm <- exp(pi * 2i / c(T1, -T1, T1, -T1))
  expectred.mu <- exp(pi * 2i / c(T2, -T2, T2, -T2))

  v1 <- sin(2 * pi * (1:N1) / T1)
  v2 <- sin(2 * pi * (1:N2) / T2)
  mx <- outer(v1, v2)
  ss <- ssa(mx, kind = "2d-ssa", L = c(N1, N2), circular = TRUE)
  for (solve.method in c("ls", "tls")) {
    for (pairing.method in c("memp", "diag")) {
      par <- parestimate(ss, groups = list(sines = 1:4),
                         solve.method = solve.method,
                         pairing.method = pairing.method)
      par <- parestimate(ss, groups = list(sines = 1:4), method = method)
      lm <- par[[1]]$roots
      mu <- par[[2]]$roots

      expect_true(is_multisets_approx_equal(mu, expectred.mu) && is_multisets_approx_equal(lm , expectred.lm),
                  label = sprintf("Est. ch. roots for sum of sines with periods %3.1f and %3.1f are correct (method = %s-%s)",
                                  T1, T2, solve.method, pairing.method))
    }
  }
})

test_that("parestimate works correctly for two sines in cylindrical 2d case", {
  r <- 4
  N1 <- 55
  N2 <- 40
  T1 <- 5
  T2 <- 8
  expectred.lm <- exp(pi * 2i / c(T1, -T1, T1, -T1))
  expectred.mu <- exp(pi * 2i / c(T2, -T2, T2, -T2))

  v1 <- sin(2 * pi * (1:N1) / T1)
  v2 <- sin(2 * pi * (1:N2) / T2)
  mx <- outer(v1, v2)
  ss <- ssa(mx, kind = "2d-ssa", L = c(N1, 20), circular = c(TRUE, FALSE))
  for (solve.method in c("ls", "tls")) {
    for (pairing.method in c("memp", "diag")) {
    par <- parestimate(ss, groups = list(sines = 1:4), method = method)
    lm <- par[[1]]$roots
    mu <- par[[2]]$roots

    expect_true(is_multisets_approx_equal(mu, expectred.mu) && is_multisets_approx_equal(lm , expectred.lm),
                label = sprintf("Est. ch. roots for sum of sines with periods %3.1f and %3.1f are correct (method = %s-%s)",
                                T1, T2, solve.method, pairing.method))
    }
  }
})
