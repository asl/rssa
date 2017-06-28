#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "extmat.h"
#include "fft_plan.h"


/* .Call calls */
extern SEXP convolveN(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hankelize_multi_fft(SEXP, SEXP, SEXP);
extern SEXP hankelize_one_fft(SEXP, SEXP, SEXP);
extern SEXP hbhankelize_one_fft(SEXP, SEXP, SEXP);
extern SEXP initialize_fft_plan(SEXP, SEXP, SEXP, SEXP);
extern SEXP initialize_hbhmat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP initialize_hmat(SEXP, SEXP, SEXP, SEXP);
extern SEXP initialize_tmat(SEXP, SEXP);
extern SEXP is_extptrnull(SEXP);
extern SEXP is_fft_plan(SEXP);
extern SEXP is_hbhmat(SEXP);
extern SEXP is_hmat(SEXP);
extern SEXP is_tmat(SEXP);
extern SEXP Lcor(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"convolveN",           (DL_FUNC) &convolveN,           5},
    {"hankelize_multi_fft", (DL_FUNC) &hankelize_multi_fft, 3},
    {"hankelize_one_fft",   (DL_FUNC) &hankelize_one_fft,   3},
    {"hbhankelize_one_fft", (DL_FUNC) &hbhankelize_one_fft, 3},
    {"initialize_fft_plan", (DL_FUNC) &initialize_fft_plan, 4},
    {"initialize_hbhmat",   (DL_FUNC) &initialize_hbhmat,   6},
    {"initialize_hmat",     (DL_FUNC) &initialize_hmat,     4},
    {"initialize_tmat",     (DL_FUNC) &initialize_tmat,     2},
    {"is_extptrnull",       (DL_FUNC) &is_extptrnull,       1},
    {"is_fft_plan",         (DL_FUNC) &is_fft_plan,         1},
    {"is_hbhmat",           (DL_FUNC) &is_hbhmat,           1},
    {"is_hmat",             (DL_FUNC) &is_hmat,             1},
    {"is_tmat",             (DL_FUNC) &is_tmat,             1},
    {"Lcor_",               (DL_FUNC) &Lcor,                3},
    {NULL, NULL, 0}
};

void R_init_Rssa(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
