#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "declarations.h"

static R_FortranMethodDef FortEntries[] = {
		{"fartransform", (DL_FUNC) &F77_SUB(artransform), 2},
		{"fldl", (DL_FUNC) &F77_SUB(ldl), 4},
		{"fldlssm", (DL_FUNC) &F77_SUB(ldlssm), 15},
		{"fsignaltheta", (DL_FUNC) &F77_SUB(signaltheta), 12},
		{"fapprox", (DL_FUNC) &F77_SUB(approx), 26},
		{"fgsmoothall", (DL_FUNC) &F77_SUB(gsmoothall), 43},
		{"fngsmooth", (DL_FUNC) &F77_SUB(ngsmooth), 39},
		{"fkfilter", (DL_FUNC) &F77_SUB(kfilter), 31},
		{"fgloglik", (DL_FUNC) &F77_SUB(gloglik), 19},
		{"fngloglik", (DL_FUNC) &F77_SUB(ngloglik), 36},
		{"fisample", (DL_FUNC) &F77_SUB(isample), 35},
		{"fzalpha", (DL_FUNC) &F77_SUB(zalpha), 10},
		{"fvarmeanw", (DL_FUNC) &F77_SUB(varmeanw), 8},
		{"fsimfilter", (DL_FUNC) &F77_SUB(simfilter), 30},
		{"fngfilter", (DL_FUNC) &F77_SUB(ngfilter), 39},
		{"fisamplefilter", (DL_FUNC) &F77_SUB(isamplefilter), 35},
		{"fsimgaussian", (DL_FUNC) &F77_SUB(simgaussian), 30},
		{"fmvfilter", (DL_FUNC) &F77_SUB(mvfilter), 12},
		{NULL, NULL, 0}
};

void R_init_KFAS(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
