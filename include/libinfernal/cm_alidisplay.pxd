from libc.stdio cimport FILE
from libc.stdint cimport int64_t, uint32_t, uint64_t

from libeasel.sq cimport ESL_SQ
from libinfernal.cm cimport CM_t


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

    # int            cm_alidisplay_Create(CM_t *cm, char *errbuf, CM_ALNDATA *adata, const ESL_SQ *sq, int64_t seqoffset, double tau, double elapsed_secs, CM_ALIDISPLAY **ret_ad)
    # int            cm_alidisplay_CreateFromP7(CM_t *cm, char *errbuf, const ESL_SQ *sq, int64_t seqoffset, float p7sc, float p7pp, P7_ALIDISPLAY *p7ad, CM_ALIDISPLAY **ret_ad)
    CM_ALIDISPLAY* cm_alidisplay_Clone(const CM_ALIDISPLAY *ad)
    size_t         cm_alidisplay_Sizeof(const CM_ALIDISPLAY *ad)
    void           cm_alidisplay_Destroy(CM_ALIDISPLAY *ad)
    char           cm_alidisplay_EncodePostProb(float p)
    float          cm_alidisplay_DecodePostProb(char pc)
    bint           cm_alidisplay_Print(FILE *fp, CM_ALIDISPLAY *ad, int min_aliwidth, int linewidth, int show_accessions)
    bint           cm_alidisplay_Is5PTrunc     (const CM_ALIDISPLAY *ad)
    bint           cm_alidisplay_Is3PTrunc     (const CM_ALIDISPLAY *ad)
    bint           cm_alidisplay_Is5PAnd3PTrunc(const CM_ALIDISPLAY *ad)
    bint           cm_alidisplay_Is5PTruncOnly (const CM_ALIDISPLAY *ad)
    bint           cm_alidisplay_Is3PTruncOnly (const CM_ALIDISPLAY *ad)
    char*          cm_alidisplay_TruncString   (const CM_ALIDISPLAY *ad)
    # int            cm_alidisplay_Backconvert(CM_t *cm, const CM_ALIDISPLAY *ad, char *errbuf, ESL_SQ **ret_sq, Parsetree_t **ret_tr, char **ret_pp)
    # int            cm_alidisplay_Dump(FILE *fp, const CM_ALIDISPLAY *ad)
    # int            cm_alidisplay_Compare(const CM_ALIDISPLAY *ad1, const CM_ALIDISPLAY *ad2)
