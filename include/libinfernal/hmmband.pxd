# from libinfernal.cm cimport CM_t

cdef extern from "infernal.h" nogil:
    
    cdef struct cp9bands_s:
        int hmm_M
        int cm_M
        int *pn_min_m
        int *pn_max_m
        int *pn_min_i
        int *pn_max_i
        int *pn_min_d
        int *pn_max_d
        int *isum_pn_m
        int *isum_pn_i
        int *isum_pn_d
        int sp1
        int ep1
        int sp2
        int ep2
        float thresh1
        float thresh2
        int Rmarg_imin
        int Rmarg_imax
        int Lmarg_jmin
        int Lmarg_jmax
        int *Jvalid
        int *Lvalid
        int *Rvalid
        int *Tvalid
        int *imin
        int *imax
        int *jmin
        int *jmax
        int **hdmin
        int **hdmax
        int *hdmin_mem
        int *hdmax_mem
        int *safe_hdmin
        int *safe_hdmax
        int hd_needed
        int hd_alloced
        double tau
    ctypedef cp9bands_s CP9Bands_t

    cdef CP9Bands_t  *AllocCP9Bands(int cm_M, int hmm_M)
    cdef float        SizeofCP9Bands(CP9Bands_t *cp9b)
    cdef void         FreeCP9Bands(CP9Bands_t *cp9bands)
    # cdef int          cp9_HMM2ijBands(CM_t *cm, char *errbuf, CP9_t *cp9, CP9Bands_t *cp9b, CP9Map_t *cp9map, int i0, int j0, int doing_search, int do_trunc, int debug_level)
    # cdef int          cp9_HMM2ijBands_OLD(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, CP9Map_t *cp9map, int i0, int j0, int doing_search, int debug_level)
    # cdef int          cp9_Seq2Bands     (CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, int doing_search, int pass_idx, int debug_level)
    # cdef int          cp9_IterateSeq2Bands(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int64_t i0, int64_t j0, int pass_idx, float size_limit, int doing_search, int do_sample, int do_post, int do_iterate, double maxtau, float *ret_Mb)
    # cdef int          cp9_Seq2Posteriors(CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int i0, int j0, int debug_level)
    # cdef void         cp9_DebugPrintHMMBands(FILE *ofp, int L, CP9Bands_t *cp9b, double hmm_bandp, int debug_level)
    # cdef int          cp9_GrowHDBands(CP9Bands_t *cp9b, char *errbuf)
    # cdef int          cp9_ValidateBands(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int i0, int j0, int do_trunc)
    # cdef void         cp9_ShiftCMBands(CM_t *cm, int i, int j, int do_trunc)
    # cdef CP9Bands_t  *cp9_CloneBands(CP9Bands_t *src_cp9b, char *errbuf)
    # cdef void         cp9_PredictStartAndEndPositions(CP9_MX *pmx, CP9Bands_t *cp9b, int i0, int j0)
    # cdef int          cp9_MarginalCandidatesFromStartEndPositions(CM_t *cm, CP9Bands_t *cp9b, int pass_idx, char *errbuf)
    # cdef void         ij2d_bands(CM_t *cm, int L, int *imin, int *imax, int *jmin, int *jmax, int **hdmin, int **hdmax, int do_trunc, int debug_level)
    # cdef void         PrintDPCellsSaved_jd(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax, int W)
    # cdef void         debug_print_ij_bands(CM_t *cm)