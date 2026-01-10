from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_profile cimport P7_PROFILE


cdef extern from "infernal.h" nogil:

    int p7_ProfileConfig5PrimeTrunc(P7_PROFILE *gm, int L)
    int p7_ProfileConfig3PrimeTrunc(const P7_HMM *hmm, P7_PROFILE *gm, int L)
    int p7_ProfileConfig5PrimeAnd3PrimeTrunc(P7_PROFILE *gm, int L)
    int p7_ReconfigLength5PrimeTrunc(P7_PROFILE *gm, int L)
    int p7_ReconfigLength3PrimeTrunc(P7_PROFILE *gm, int L)
