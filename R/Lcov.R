#Useless function. Will be deleted soon
#We are not needed to get fist row of Lcov matrix. 
#Lcov <- function(F, L) {
#  storage.mode(F) <- "double";
#  storage.mode(L) <- "integer";
#  .Call("Lcov", F, L);
#}

Lcov.matrix <- function(F, L) {
  storage.mode(F) <- "double";
  storage.mode(L) <- "integer";
  .Call("Lcov_matrix", F, L);
}
