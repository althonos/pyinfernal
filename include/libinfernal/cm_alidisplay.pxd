from libc.stdint cimport int64_t, uint32_t, uint64_t


cdef extern from "infernal.h" nogil:

    cdef struct cm_alidisplay_s:
        char *rfline
        char *ncline
        char *csline
        char *model
        char *mline
        char *aseq
        char *ppline
        int   N
        char *aseq_el
        char *rfline_el
        char *ppline_el
        int   N_el
        char *cmname
        char *cmacc
        char *cmdesc
        int   cfrom_emit
        int   cto_emit
        int   cfrom_span
        int   cto_span
        int   clen
        char *sqname
        char *sqacc
        char *sqdesc
        long  sqfrom
        long  sqto
        float  sc
        float  avgpp
        float  gc
        double tau
        float  matrix_Mb
        double elapsed_secs
        bint   hmmonly
        int   memsize
        char *mem
    ctypedef cm_alidisplay_s CM_ALIDISPLAY