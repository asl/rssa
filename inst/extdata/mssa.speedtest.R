library(Rssa)
library(microbenchmark)


mssa <- function(v, L, r,
                 method = c("2d-ssa", "mssa"),
                 svd.method = "nutrlan") {
  method <- match.arg(method)

  N <- nrow(v)
  K <- (N - L + 1) * ncol(v)
  rr <- min(L, K, r + 1)

  if (identical(method, "2d-ssa"))
    ssa(v, L = c(L, 1), kind = "2d-ssa", svd.method = svd.method, neig = rr)
  else if(identical(method, "mssa"))
    ssa(v, L = L, kind = "mssa", svd.method = svd.method, neig = rr)
  else
    stop("Unknown `method'")
}

err_mssa <- function(v, L, r, sigma, ...) {
  N <- nrow(v)
  s <- mssa(v + rnorm(prod(dim(v)), sd = sigma), L, r, ...)
  H <- reconstruct(s, groups = list(1:r), drop.attributes = TRUE)$F1
  mean((as.numeric(v) - as.numeric(unlist(H)))^2, na.rm = TRUE)
}

test <- function(sigma = 1,
                 Nser = 119,
                 r = 2,
                 nums = 20,
                 R = 20,
                 verbose = FALSE,
                 ...) {
  dots <- list(...)
  res <- sapply(nums, function(num) {
      if (verbose) cat(sprintf("num = %d\n", num))

      v <- sapply(seq_len(num), function(i) cos(seq_len(Nser) * 2*pi/10 + i*pi/num))
      Ls <- seq(10, 110, 10)
      res <- sapply(Ls, function(L) {
          if (verbose) cat(sprintf("L = %d; ", L))
          err <- replicate(R, do.call("err_mssa", c(list(v, L, r, sigma), dots)))
          mean(err)
        })
      names(res) <- as.character(Ls)

      if (verbose) cat("\n")
      res
    })

  res
}

cat("MSSA:\n")
set.seed(1)
print(microbenchmark(test(method = "mssa"), times = 5))

cat("2d-SSA:\n")
set.seed(1)
print(microbenchmark(test(method = "2d-ssa"), times = 5))
