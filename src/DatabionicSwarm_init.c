#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _DatabionicSwarm_Delta3DWeightsC(SEXP, SEXP);
extern SEXP _DatabionicSwarm_DijkstraSSSP(SEXP, SEXP, SEXP);
extern SEXP _DatabionicSwarm_findPossiblePositionsCsingle(SEXP, SEXP, SEXP, SEXP);
extern SEXP _DatabionicSwarm_PswarmCurrentRadiusC2botsPositive(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DatabionicSwarm_rDistanceToroidCsingle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DatabionicSwarm_trainstepC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DatabionicSwarm_trainstepC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"_DatabionicSwarm_Delta3DWeightsC",                   (DL_FUNC) &_DatabionicSwarm_Delta3DWeightsC,                    2},
    {"_DatabionicSwarm_DijkstraSSSP",                      (DL_FUNC) &_DatabionicSwarm_DijkstraSSSP,                       3},
    {"_DatabionicSwarm_findPossiblePositionsCsingle",      (DL_FUNC) &_DatabionicSwarm_findPossiblePositionsCsingle,       4},
    {"_DatabionicSwarm_PswarmCurrentRadiusC2botsPositive", (DL_FUNC) &_DatabionicSwarm_PswarmCurrentRadiusC2botsPositive, 14},
    {"_DatabionicSwarm_rDistanceToroidCsingle",            (DL_FUNC) &_DatabionicSwarm_rDistanceToroidCsingle,             6},
    {"_DatabionicSwarm_trainstepC",                        (DL_FUNC) &_DatabionicSwarm_trainstepC,                         8},
    {"_DatabionicSwarm_trainstepC2",                       (DL_FUNC) &_DatabionicSwarm_trainstepC2,                        9},
    {NULL, NULL, 0}
};

void R_init_DatabionicSwarm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}