from libc.stdint cimport int64_t, uint32_t, uint64_t
from libc.stdio cimport FILE

from libeasel.keyhash cimport ESL_KEYHASH
from libinfernal.cm cimport CM_t
from libinfernal.cm_alidisplay cimport CM_ALIDISPLAY


cdef extern from "infernal.h" nogil:

    cdef struct cm_pipeline_s:
        pass
    ctypedef cm_pipeline_s CM_PIPELINE

    cdef struct cm_hit_s:
        char          *name
        char          *acc
        char          *desc
        int64_t        cm_idx
        int            clan_idx
        int64_t        seq_idx
        int            pass_idx
        int64_t        hit_idx
        int64_t        srcL
        int64_t        start
        int64_t        stop
        bint           in_rc
        int            root
        int            mode
        float          score
        float          bias
        double         pvalue
        double         evalue
        bint           has_evalue
        bint           hmmonly
        bint           glocal
        CM_ALIDISPLAY *ad
        uint32_t       flags
        int64_t        any_oidx
        int64_t        win_oidx
        double         any_bitE
        double         win_bitE
    ctypedef cm_hit_s CM_HIT

    cdef struct cm_tophits_s:
        CM_HIT **hit         
        CM_HIT  *unsrt	      
        uint64_t Nalloc
        uint64_t N
        uint64_t nreported
        uint64_t nincluded
        bint      is_sorted_by_evalue
        bint      is_sorted_for_overlap_removal
        bint      is_sorted_for_overlap_markup
        bint      is_sorted_by_position
    ctypedef cm_tophits_s CM_TOPHITS

    cdef enum:
        CM_HIT_FLAGS_DEFAULT 
        CM_HIT_IS_INCLUDED            
        CM_HIT_IS_REPORTED            
        CM_HIT_IS_REMOVED_DUPLICATE   
        CM_HIT_IS_MARKED_OVERLAP      

    cdef CM_TOPHITS *cm_tophits_Create()
    cdef int         cm_tophits_Grow(CM_TOPHITS *h)
    cdef int         cm_tophits_CreateNextHit(CM_TOPHITS *h, CM_HIT **ret_hit)
    cdef int         cm_tophits_SortByEvalue(CM_TOPHITS *h)
    cdef int         cm_tophits_SortForOverlapRemoval(CM_TOPHITS *h)
    cdef int         cm_tophits_SortForOverlapMarkup(CM_TOPHITS *h, int do_clans_only)
    cdef int         cm_tophits_SortByPosition(CM_TOPHITS *h)
    cdef int         cm_tophits_Merge(CM_TOPHITS *h1, CM_TOPHITS *h2)
    cdef int         cm_tophits_GetMaxPositionLength(CM_TOPHITS *h)
    cdef int         cm_tophits_GetMaxTargetLength(CM_TOPHITS *h)
    cdef int         cm_tophits_GetMaxNameLength(CM_TOPHITS *h)
    cdef int         cm_tophits_GetMaxDescLength(CM_TOPHITS *h)
    cdef int         cm_tophits_GetMaxAccessionLength(CM_TOPHITS *h)
    cdef int         cm_tophits_GetMaxShownLength(CM_TOPHITS *h)
    cdef int         cm_tophits_GetMaxClanLength(CM_TOPHITS *h, ESL_KEYHASH *clan_name_kh)
    cdef int         cm_tophits_GetMaxModelLength(CM_TOPHITS *h)
    cdef int         cm_tophits_Reuse(CM_TOPHITS *h)
    cdef void        cm_tophits_Destroy(CM_TOPHITS *h)
    cdef int         cm_tophits_CloneHitMostly(CM_TOPHITS *src_th, int h, CM_TOPHITS *dest_th)
    cdef int         cm_tophits_ComputeEvalues(CM_TOPHITS *th, double eZ, int istart)
    cdef int         cm_tophits_RemoveOrMarkOverlaps(CM_TOPHITS *th, int do_clans_only, char *errbuf)
    cdef int         cm_tophits_UpdateHitPositions(CM_TOPHITS *th, int hit_start, int64_t seq_start, bint in_revcomp)
    cdef int         cm_tophits_SetSourceLengths(CM_TOPHITS *th, int64_t *srcL, uint64_t nseqs)
    cdef int64_t     cm_tophits_OverlapNres(int64_t from1, int64_t to1, int64_t from2, int64_t to2, int64_t *ret_nes, char *errbuf)

    cdef int cm_tophits_Threshold(CM_TOPHITS *th, CM_PIPELINE *pli)
    # cdef int cm_tophits_Targets(FILE *ofp, CM_TOPHITS *th, CM_PIPELINE *pli, int textw)
    # cdef int cm_tophits_F3Targets(FILE *ofp, CM_TOPHITS *th, CM_PIPELINE *pli)
    # cdef int cm_tophits_HitAlignments(FILE *ofp, CM_TOPHITS *th, CM_PIPELINE *pli, int textw)
    # cdef int cm_tophits_HitAlignmentStatistics(FILE *ofp, CM_TOPHITS *th, int used_cyk, int used_hb, double default_tau)
    # cdef int cm_tophits_Alignment(CM_t *cm, const CM_TOPHITS *th, char *errbuf, int allow_trunc, ESL_MSA **ret_msa)
    cdef int cm_tophits_TabularTargets1(FILE *ofp, char *qname, char *qacc, CM_TOPHITS *th, CM_PIPELINE *pli, int show_header) except *
    cdef int cm_tophits_TabularTargets2(FILE *ofp, char *qname, char *qacc, CM_TOPHITS *th, CM_PIPELINE *pli, int show_header, ESL_KEYHASH *clan_name_kh, int skip_overlaps, char *errbuf) except *
    cdef int cm_tophits_TabularTargets3(FILE *ofp, char *qname, char *qacc, CM_TOPHITS *th, CM_PIPELINE *pli, int show_header) except *
    # cdef int cm_tophits_F3TabularTargets1(FILE *ofp, CM_TOPHITS *th, CM_PIPELINE *pli, int show_header)
    # cdef int cm_tophits_TabularTail(FILE *ofp, const char *progname, enum cm_pipemodes_e pipemode, const char *qfile, const char *tfile, const ESL_GETOPTS *go)
    # cdef int cm_tophits_Dump(FILE *fp, const CM_TOPHITS *th)

    cdef int    cm_hit_AllowTruncation(CM_t *cm, int pass_idx, int64_t start, int64_t stop, int64_t i0, int64_t j0, char mode, int b)
    cdef int    cm_hit_Dump(FILE *fp, const CM_HIT *h)
