from libinfernal.cm cimport CM_t


cdef extern from "infernal.h" nogil:

    int   cm_Configure(CM_t *cm, char *errbuf, int W_from_cmdline)
    # int   cm_ConfigureSub(CM_t *cm, char *errbuf, int W_from_cmdline, CM_t *mother_cm, CMSubMap_t *mother_map)
    # int   cm_CalculateLocalBeginProbs(CM_t *cm, float p_internal_start, float **t, float *begin)
