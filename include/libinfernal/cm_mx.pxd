from libc.stdint cimport int64_t
from libc.stdio cimport FILE

from libinfernal.cm_qdband cimport CM_QDBINFO
from libinfernal.hmmband cimport CP9Bands_t


cdef extern from "infernal.h" nogil:

    cdef struct cm_s:
        pass
    ctypedef cm_s CM_t
        
    cdef struct cm_mx_s:
        int M
        int L
        int64_t ncells_alloc
        int64_t ncells_valid
        float size_Mb
        float*** dp
        float* dp_mem
    ctypedef cm_mx_s CM_MX

    cdef struct cm_tr_mx_s:
        int M
        int B
        int L
        int64_t Jncells_alloc
        int64_t Jncells_valid
        int64_t Lncells_alloc
        int64_t Lncells_valid
        int64_t Rncells_alloc
        int64_t Rncells_valid
        int64_t Tncells_alloc
        int64_t Tncells_valid
        float size_Mb
        float*** Jdp
        float* Jdp_mem
        float*** Ldp
        float* Ldp_mem
        float*** Rdp
        float* Rdp_mem
        float*** Tdp
        float* Tdp_mem
    ctypedef cm_tr_mx_s CM_TR_MX

    cdef struct cm_hb_mx_s:
        int  M
        int  L
        int64_t ncells_alloc
        int64_t ncells_valid
        float size_Mb
        int* nrowsA
        float*** dp
        float* dp_mem
        CP9Bands_t* cp9b
    ctypedef cm_hb_mx_s CM_HB_MX

    cdef struct cm_tr_hb_mx_s:
        int M
        int B
        int L
        int64_t Jncells_alloc
        int64_t Lncells_alloc
        int64_t Rncells_alloc
        int64_t Tncells_alloc
        int64_t Jncells_valid
        int64_t Lncells_valid
        int64_t Rncells_valid
        int64_t Tncells_valid
        float size_Mb
        int* JnrowsA
        int* LnrowsA
        int* RnrowsA
        int* TnrowsA
        float*** Jdp
        float* Jdp_mem
        float*** Ldp
        float* Ldp_mem
        float*** Rdp
        float* Rdp_mem
        float*** Tdp
        float* Tdp_mem
        CP9Bands_t* cp9b
    cdef cm_tr_hb_mx_s CM_TR_HB_MX

    cdef struct cm_shadow_mx_s:
        int M
        int L
        int B
        # int64_t y_ncells_alloc;  /* current cell allocation limit in yshadow*/
        # int64_t k_ncells_alloc;  /* current cell allocation limit in kshadow*/
        # int64_t y_ncells_valid;  /* current number of valid cells in yshadow */
        # int64_t k_ncells_valid;  /* current number of valid cells in kshadow */
        # float  size_Mb;         /* current size of matrix in Megabytes */
        # char ***yshadow;       /* [0..v..M][0..j..L][0..d..j] */
        # char   *yshadow_mem;   /* the actual mem, points to yshadow[0][0][0] */
        # int ***kshadow;       /* [0..v..M][0..j..L][0..d..j] */
        # int   *kshadow_mem;   /* the actual mem, points to Jkshadow[0][0][0] */
    ctypedef cm_shadow_mx_s CM_SHADOW_MX

    cdef struct cm_tr_shadow_mx_s:
        int  M
        int  L
        int  B
        # int64_t Jy_ncells_alloc;  /* current cell allocation limit in Jyshadow*/
        # int64_t Ly_ncells_alloc;  /* current cell allocation limit in Lyshadow*/
        # int64_t Ry_ncells_alloc;  /* current cell allocation limit in Ryshadow*/
        # int64_t Jk_ncells_alloc;  /* current cell allocation limit in Jkshadow*/
        # int64_t Lk_ncells_alloc;  /* current cell allocation limit in Lkshadow/Lkmode*/
        # int64_t Rk_ncells_alloc;  /* current cell allocation limit in Rkshadow/Rkmode*/
        # int64_t Tk_ncells_alloc;  /* current cell allocation limit in Tkshadow*/
        # int64_t Jy_ncells_valid;  /* current number of valid cells in Jyshadow */
        # int64_t Ly_ncells_valid;  /* current number of valid cells in Lyshadow */
        # int64_t Ry_ncells_valid;  /* current number of valid cells in Ryshadow */
        # int64_t Jk_ncells_valid;  /* current number of valid cells in Jkshadow */
        # int64_t Lk_ncells_valid;  /* current number of valid cells in Lkshadow/Lkmode */
        # int64_t Rk_ncells_valid;  /* current number of valid cells in Rkshadow/Rkmode */
        # int64_t Tk_ncells_valid;  /* current number of valid cells in Tkshadow */
        # float  size_Mb;         /* current size of matrix in Megabytes */
        # char ***Jyshadow;       /* [0..v..M][0..j..L][0..d..j] */
        # char   *Jyshadow_mem;   /* the actual mem, points to Jyshadow[0][0][0] */
        # char ***Lyshadow;       /* [0..v..M][0..j..L][0..d..j] */
        # char   *Lyshadow_mem;   /* the actual mem, points to Lyshadow[0][0][0] */
        # char ***Ryshadow;       /* [0..v..M][0..j..L][0..d..j] */
        # char   *Ryshadow_mem;   /* the actual mem, points to Ryshadow[0][0][0] */
        # int ***Jkshadow;       /* [0..v..M][0..j..L][0..d..j] */
        # int   *Jkshadow_mem;   /* the actual mem, points to Jkshadow[0][0][0] */
        # int ***Lkshadow;       /* [0..v..M][0..j..L][0..d..j] */
        # int   *Lkshadow_mem;   /* the actual mem, points to Lkshadow[0][0][0] */
        # int ***Rkshadow;       /* [0..v..M][0..j..L][0..d..j] */
        # int   *Rkshadow_mem;   /* the actual mem, points to Rkshadow[0][0][0] */
        # int ***Tkshadow;       /* [0..v..M][0..j..L][0..d..j] */
        # int   *Tkshadow_mem;   /* the actual mem, points to Tkshadow[0][0][0] */
        # char ***Lkmode;        /* [0..v..M][0..j..L][0..d..j] */
        # char   *Lkmode_mem;    /* the actual mem, points to Lkmode[0][0][0] */
        # char ***Rkmode;        /* [0..v..M][0..j..L][0..d..j] */
        # char   *Rkmode_mem;    /* the actual mem, points to Rkmode[0][0][0] */
    ctypedef cm_tr_shadow_mx_s CM_TR_SHADOW_MX

    cdef struct cm_hb_shadow_mx_s:
        int M
        int L
        int B
        int64_t y_ncells_alloc
        int64_t y_ncells_valid
        int64_t k_ncells_alloc
        int64_t k_ncells_valid
        float size_Mb
        int* nrowsA
        char*** yshadow
        char* yshadow_mem
        int*** kshadow
        int* kshadow_mem
        CP9Bands_t *cp9b
    ctypedef cm_hb_shadow_mx_s CM_HB_SHADOW_MX

    cdef struct cm_tr_hb_shadow_mx_s:
        int M
        int L
        int B
        # int64_t Jy_ncells_alloc;  /* current cell allocation limit in Jyshadow*/
        # int64_t Ly_ncells_alloc;  /* current cell allocation limit in Lyshadow*/
        # int64_t Ry_ncells_alloc;  /* current cell allocation limit in Ryshadow*/
        # int64_t Jk_ncells_alloc;  /* current cell allocation limit in Jkshadow*/
        # int64_t Lk_ncells_alloc;  /* current cell allocation limit in Lkshadow/Lkmode*/
        # int64_t Rk_ncells_alloc;  /* current cell allocation limit in Rkshadow/Rkmode*/
        # int64_t Tk_ncells_alloc;  /* current cell allocation limit in Tkshadow*/
        # int64_t Jy_ncells_valid;  /* current number of valid cells in Jyshadow */
        # int64_t Ly_ncells_valid;  /* current number of valid cells in Lyshadow */
        # int64_t Ry_ncells_valid;  /* current number of valid cells in Ryshadow */
        # int64_t Jk_ncells_valid;  /* current number of valid cells in Jkshadow */
        # int64_t Lk_ncells_valid;  /* current number of valid cells in Lkshadow/Lkmode */
        # int64_t Rk_ncells_valid;  /* current number of valid cells in Rkshadow/Rkmode */
        # int64_t Tk_ncells_valid;  /* current number of valid cells in Tkshadow */
        # float  size_Mb;         /* current size of matrix in Megabytes */
        # int   *JnrowsA;         /* [0..v..M] current number allocated rows for deck v */
        # int   *LnrowsA;         /* [0..v..M] current number allocated rows for deck v */
        # int   *RnrowsA;         /* [0..v..M] current number allocated rows for deck v */
        # int   *TnrowsA;         /* [0..v..M] current number allocated rows for deck v */
        # char ***Jyshadow;       /* [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
        # char   *Jyshadow_mem;   /* the actual mem, points to Jyshadow[0][0][0] */
        # char ***Lyshadow;       /* [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
        # char   *Lyshadow_mem;   /* the actual mem, points to Lyshadow[0][0][0] */
        # char ***Ryshadow;       /* [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
        # char   *Ryshadow_mem;   /* the actual mem, points to Ryshadow[0][0][0] */
        # int ***Jkshadow;       /*  [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
        # int   *Jkshadow_mem;   /* the actual mem, points to Jkshadow[0][0][0] */
        # int ***Lkshadow;       /*  [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
        # int   *Lkshadow_mem;   /* the actual mem, points to Lkshadow[0][0][0] */
        # int ***Rkshadow;       /*  [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
        # int   *Rkshadow_mem;   /* the actual mem, points to Rkshadow[0][0][0] */
        # int ***Tkshadow;       /*  [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
        # int   *Tkshadow_mem;   /* the actual mem, points to Tkshadow[0][0][0] */
        # char ***Lkmode;        /*  [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
        # char   *Lkmode_mem;    /* the actual mem, points to Lkmode[0][0][0] */
        # char ***Rkmode;        /*  [0..v..M][0..j..(cp9b->jmax[v]-cp9b->jmin[v])[0..d..cp9b->hdmax[v][j-jmin[v]]-cp9b->hdmin[v][j-jmin[v]]] */
        # char   *Rkmode_mem;    /* the actual mem, points to Rkmode[0][0][0] */
        # CP9Bands_t *cp9b;     
    ctypedef cm_tr_hb_shadow_mx_s CM_TR_HB_SHADOW_MX

    cdef struct cm_emit_mx_s:
        int M
        int L
        # int64_t  l_ncells_alloc;  /* current cell allocation limit for l_pp */
        # int64_t  l_ncells_valid;  /* current number of valid cells for l_pp */
        # int64_t  r_ncells_alloc;  /* current cell allocation limit for r_pp */
        # int64_t  r_ncells_valid;  /* current number of valid cells for r_pp */
        # float    size_Mb;         /* current size of matrix in Megabytes  */
        # float   **l_pp;         /* matrix: [0..v..M][0..1..i..L], l_pp[v][0] is
        #             * always IMPOSSIBLE l_pp[v] == NULL if v is
        #             * not a left emitter.
        #             */
        # float   **r_pp;         /* matrix: [0..v..M][0..1..i..L], r_pp[v][0] is
        #             * always IMPOSSIBLE r_pp[v] == NULL if v is
        #             * not a right emitter.
        #             */
        # float    *l_pp_mem;     /* the actual mem for l_pp, points to
        #             * l_pp[v][0], where v is min v for which
        #             * l_pp != NULL */
        # float    *r_pp_mem;     /* the actual mem for r_pp, points to
        #             * r_pp[v][0], where v is min v for which
        #             * r_pp != NULL */
        # float    *sum;          /* [0..1..i..L] log of the summed posterior
        #             * probability that residue i was emitted
        #             * either leftwise or rightwise by any state.
        #             * Used for normalizing l_pp and r_pp.
        #             */
    ctypedef cm_emit_mx_s CM_EMIT_MX

    cdef struct cm_tr_emit_mx_s:
        int M
        int L
        # int64_t  l_ncells_alloc;   /* current cell allocation limit for Jl_pp, Ll_pp */
        # int64_t  l_ncells_valid;   /* current number of valid cells for Jl_pp, Ll_pp */
        # int64_t  r_ncells_alloc;   /* current cell allocation limit for Jr_pp, Rr_pp */
        # int64_t  r_ncells_valid;   /* current number of valid cells for Jr_pp, Rr_pp */
        # float    size_Mb;          /* current size of matrix in Megabytes  */
        # float   **Jl_pp;         /* matrix: [0..v..M][0..1..i..L], Joint mode */
        # float   **Ll_pp;         /* matrix: [0..v..M][0..1..i..L], Left mode */
        # float   **Jr_pp;         /* matrix: [0..v..M][0..1..i..L], Joint mode */
        # float   **Rr_pp;         /* matrix: [0..v..M][0..1..i..L], Right mode */
        # float    *Jl_pp_mem;     /* the actual mem for Jl_pp */
        # float    *Ll_pp_mem;     /* the actual mem for Ll_pp */
        # float    *Jr_pp_mem;     /* the actual mem for Jr_pp */
        # float    *Rr_pp_mem;     /* the actual mem for Rr_pp */
        # float    *sum;  
    ctypedef cm_tr_emit_mx_s CM_TR_EMIT_MX

    cdef struct cm_hb_emit_mx_s:
        int M
        int L
        # int64_t  l_ncells_alloc;  /* current cell allocation limit for dp */
        # int64_t  l_ncells_valid;  /* current number of valid cells for dp */
        # int64_t  r_ncells_alloc;  /* current cell allocation limit for dp */
        # int64_t  r_ncells_valid;  /* current number of valid cells for dp */
        # float    size_Mb;         /* current size of matrix in Megabytes  */
        # float   **l_pp
        # float   **r_pp
        # float    *l_pp_mem
        # float    *r_pp_mem
        # float    *sum
        # CP9Bands_t *cp9b
    ctypedef cm_hb_emit_mx_s CM_HB_EMIT_MX

    cdef struct cm_tr_hb_emit_mx_s:
        int M
        int L
        # int64_t  l_ncells_alloc;   /* current cell allocation limit for Jl_pp, Ll_pp */
        # int64_t  l_ncells_valid;   /* current number of valid cells for Jl_pp, Ll_pp */
        # int64_t  r_ncells_alloc;   /* current cell allocation limit for Jr_pp, Rr_pp */
        # int64_t  r_ncells_valid;   /* current number of valid cells for Jr_pp, Rr_pp */
        # float    size_Mb;          /* current size of matrix in Megabytes  */
        # float   **Jl_pp;         /* matrix: [0..v..M][0..1..i..L], Joint mode */
        # float   **Ll_pp;         /* matrix: [0..v..M][0..1..i..L], Left mode */
        # float   **Jr_pp;         /* matrix: [0..v..M][0..1..i..L], Joint mode */
        # float   **Rr_pp;         /* matrix: [0..v..M][0..1..i..L], Right mode */
        # float    *Jl_pp_mem;     /* the actual mem for Jl_pp */
        # float    *Ll_pp_mem;     /* the actual mem for Ll_pp */
        # float    *Jr_pp_mem;     /* the actual mem for Jr_pp */
        # float    *Rr_pp_mem;     /* the actual mem for Rr_pp */
        # float    *sum
        # CP9Bands_t *cp9b
    ctypedef cm_tr_hb_emit_mx_s CM_TR_HB_EMIT_MX

    const size_t SMX_NOQDB
    const size_t SMX_QDB1_TIGHT
    const size_t SMX_QDB2_LOOSE
    const size_t NSMX_QDB_IDX

    cdef struct cm_scan_mx_s:
        CM_QDBINFO* qdbinfo
        int M
        int W
        # int   ***dnAAA;        /* [0..n..NSMX_QDB_IDX-1][1..j..W][0..v..M-1] min d value allowed for posn j, state v using QDB set n */
        # int   ***dxAAA;        /* [0..n..NSMX_QDB_IDX-1][1..j..W][0..v..M-1] max d value allowed for posn j, state v using QDB set n */
        # int     *bestr;        /* auxil info: best root state v at alpha[0][cur][d] (0->v local begin used if v != 0)*/
        # float   *bestsc;       /* auxil info: best score for parsetree at alpha[0][cur][d] in mode bestmode[d] */
        # int      floats_valid; /* TRUE if float alpha matrices are valid, FALSE if not */
        # int      ints_valid;   /* TRUE if int   alpha matrices are valid, FALSE if not */
        # float    size_Mb;      /* size of matrix in Megabytes */
        # float ***falpha;          /* non-BEGL_S states for float versions of CYK/Inside */
        # float ***falpha_begl;     /*     BEGL_S states for float versions of CYK/Inside */
        # float   *falpha_mem;      /* ptr to the actual memory for falpha */
        # float   *falpha_begl_mem; /* ptr to the actual memory for falpha_begl */
        # int   ***ialpha;          /* non-BEGL_S states for int   versions of CYK/Inside */
        # int   ***ialpha_begl;     /*     BEGL_S states for int   versions of CYK/Inside */
        # int     *ialpha_mem;      /* ptr to the actual memory for ialpha */
        # int     *ialpha_begl_mem; /* ptr to the actual memory for ialpha_begl */
        # int64_t  ncells_alpha;      /* number of alloc'ed, valid cells for falpha and ialpha matrices, alloc'ed as contiguous block */
        # int64_t  ncells_alpha_begl; /* number of alloc'ed, valid cells for falpha_begl and ialpha_begl matrices, alloc'ed as contiguous block */
    ctypedef cm_scan_mx_s CM_SCAN_MX

    cdef struct cm_tr_scan_mx_s:
        CM_QDBINFO* qdbinfo
        int M
        int W
        # int   ***dnAAA;        /* [0..n..NSMX_QDB_IDX-1][1..j..W][0..v..M-1] min d value allowed for posn j, state v using QDB set n */
        # int   ***dxAAA;        /* [0..n..NSMX_QDB_IDX-1][1..j..W][0..v..M-1] max d value allowed for posn j, state v using QDB set n */
        # int     *bestr;        /* auxil info: best root state v at alpha[0][cur][d] (0->v truncated begin used if v != 0)*/
        # float   *bestsc;       /* auxil info: best score for parsetree at alpha[0][cur][d] in mode bestmode[d] */
        # char    *bestmode;     /* auxil info: best mode for parsetree at alpha[0][cur][d], gives score in bestsc[d] */
        # int      floats_valid; /* TRUE if float alpha matrices are valid, FALSE if not */
        # int      ints_valid;   /* TRUE if int   alpha matrices are valid, FALSE if not */
        # float    size_Mb;      /* size of matrix in Megabytes */
        # float ***fJalpha;          /* non-BEGL_S states for float versions of CYK/Inside */
        # float ***fJalpha_begl;     /*     BEGL_S states for float versions of CYK/Inside */
        # float   *fJalpha_mem;      /* ptr to the actual memory for fJalpha */
        # float   *fJalpha_begl_mem; /* ptr to the actual memory for fJalpha_begl */
        # float ***fLalpha;          /* non-BEGL_S states for float versions of CYK/Inside */
        # float ***fLalpha_begl;     /*     BEGL_S states for float versions of CYK/Inside */
        # float   *fLalpha_mem;      /* ptr to the actual memory for fLalpha */
        # float   *fLalpha_begl_mem; /* ptr to the actual memory for fLalpha_begl */
        # float ***fRalpha;          /* non-BEGL_S states for float versions of CYK/Inside */
        # float ***fRalpha_begl;     /*     BEGL_S states for float versions of CYK/Inside */
        # float   *fRalpha_mem;      /* ptr to the actual memory for fRalpha */
        # float   *fRalpha_begl_mem; /* ptr to the actual memory for fRalpha_begl */
        # float ***fTalpha;          /* BIF states for float versions of CYK/Inside */
        # float   *fTalpha_mem;      /* ptr to the actual memory for fTalpha */
        # int   ***iJalpha;          /* non-BEGL_S states for int   versions of CYK/Inside */
        # int   ***iJalpha_begl;     /*     BEGL_S states for int   versions of CYK/Inside */
        # int     *iJalpha_mem;      /* ptr to the actual memory for iJalpha */
        # int     *iJalpha_begl_mem; /* ptr to the actual memory for iJalpha_begl */
        # int   ***iLalpha;          /* non-BEGL_S states for int   versions of CYK/Inside */
        # int   ***iLalpha_begl;     /*     BEGL_S states for int   versions of CYK/Inside */
        # int     *iLalpha_mem;      /* ptr to the actual memory for iLalpha */
        # int     *iLalpha_begl_mem; /* ptr to the actual memory for iLalpha_begl */
        # int   ***iRalpha;          /* non-BEGL_S states for int   versions of CYK/Inside */
        # int   ***iRalpha_begl;     /*     BEGL_S states for int   versions of CYK/Inside */
        # int     *iRalpha_mem;      /* ptr to the actual memory for iRalpha */
        # int     *iRalpha_begl_mem; /* ptr to the actual memory for iRalpha_begl */
        # int   ***iTalpha;          /* BIF states for int   versions of CYK/Inside */
        # int     *iTalpha_mem;      /* ptr to the actual memory for iTalpha */
        # int64_t  ncells_alpha;      /* number of alloc'ed, valid cells for f{J,L,R}alpha and i{J,L,R}alpha matrices, alloc'ed as contiguous block */
        # int64_t  ncells_alpha_begl; /* number of alloc'ed, valid cells for f{J,L,R}alpha_begl and i{J,L,R}alpha_begl matrices, alloc'ed as contiguous block */
        # int64_t  ncells_Talpha;     /* number of alloc'ed, valid cells for fTalpha and iTalpha matrices, alloc'ed as contiguous block */
    ctypedef cm_tr_scan_mx_s CM_TR_SCAN_MX

    const int TRPENALTY_5P_AND_3P
    const int TRPENALTY_5P_ONLY
    const int TRPENALTY_3P_ONLY
    const int NTRPENALTY
    const int TRPENALTY_NONE

    cdef struct cm_tr_penalties_s:
        int M
        bint ignored_inserts
        float **g_ptyAA
        float **l_ptyAA
        int  **ig_ptyAA
        int  **il_ptyAA
    ctypedef cm_tr_penalties_s CM_TR_PENALTIES

    # cdef CM_MX           *cm_mx_Create                  (int M)
    # cdef int              cm_mx_GrowTo                  (CM_t *cm, CM_MX *mx, char *errbuf, int L, float size_limit)
    # cdef int              cm_mx_Dump                    (FILE *ofp, CM_MX *mx, int print_mx)
    # cdef void             cm_mx_Destroy                 (CM_MX *mx)
    # cdef int              cm_mx_SizeNeeded              (CM_t *cm, char *errbuf, int L, int64_t *ret_ncells, float *ret_Mb)

    # cdef CM_TR_MX        *cm_tr_mx_Create               (CM_t *cm)
    # cdef int              cm_tr_mx_GrowTo               (CM_t *cm, CM_TR_MX *mx, char *errbuf, int L, float size_limit)
    # cdef int              cm_tr_mx_Dump                 (FILE *ofp, CM_TR_MX *mx, char mode, int print_mx)
    # cdef void             cm_tr_mx_Destroy              (CM_TR_MX *mx)
    # cdef int              cm_tr_mx_SizeNeeded           (CM_t *cm, char *errbuf, int L, int64_t *ret_Jncells, int64_t *ret_Lncells, int64_t *ret_Rncells, int64_t *ret_Tncells, float *ret_Mb)

    # cdef CM_HB_MX        *cm_hb_mx_Create               (int M)
    # cdef int              cm_hb_mx_GrowTo               (CM_t *cm, CM_HB_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit)
    # cdef int              cm_hb_mx_Dump                 (FILE *ofp, CM_HB_MX *mx, int print_mx)
    # cdef void             cm_hb_mx_Destroy              (CM_HB_MX *mx)
    # cdef int              cm_hb_mx_SizeNeeded           (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int L, int64_t *ret_ncells, float *ret_Mb)

    # cdef CM_TR_HB_MX     *cm_tr_hb_mx_Create            (CM_t *cm)
    # cdef int              cm_tr_hb_mx_GrowTo            (CM_t *cm, CM_TR_HB_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit)
    # cdef int              cm_tr_hb_mx_Dump              (FILE *ofp, CM_TR_HB_MX *mx, char mode, int print_mx)
    # cdef void             cm_tr_hb_mx_Destroy           (CM_TR_HB_MX *mx)
    # cdef int              cm_tr_hb_mx_SizeNeeded        (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int L, int64_t *ret_Jncells, int64_t *ret_Lncells, 
    #                             int64_t *ret_Rncells, int64_t *ret_Tncells, float *ret_Mb)

    # cdef CM_SHADOW_MX    *cm_shadow_mx_Create           (CM_t *cm)
    # cdef int              cm_shadow_mx_GrowTo           (CM_t *cm, CM_SHADOW_MX *mx, char *errbuf, int L, float size_limit)
    # cdef int              cm_shadow_mx_Dump             (FILE *ofp, CM_t *cm, CM_SHADOW_MX *mx, int print_mx)
    # cdef void             cm_shadow_mx_Destroy          (CM_SHADOW_MX *mx)
    # cdef int              cm_shadow_mx_SizeNeeded       (CM_t *cm, char *errbuf, int L, int64_t *ret_ny_cells, int64_t *ret_nk_cells, float *ret_Mb)

    # cdef CM_TR_SHADOW_MX *cm_tr_shadow_mx_Create        (CM_t *cm)
    # cdef int              cm_tr_shadow_mx_GrowTo        (CM_t *cm, CM_TR_SHADOW_MX *mx, char *errbuf, int L, float size_limit)
    # cdef int              cm_tr_shadow_mx_Dump          (FILE *ofp, CM_t *cm, CM_TR_SHADOW_MX *mx, char mode, int print_mx)
    # cdef void             cm_tr_shadow_mx_Destroy       (CM_TR_SHADOW_MX *mx)
    # cdef int              cm_tr_shadow_mx_SizeNeeded    (CM_t *cm, char *errbuf, int L, int64_t *ret_Jny_cells, int64_t *ret_Lny_cells, int64_t *ret_Rny_cells,
    #                             int64_t *ret_Jnk_cells, int64_t *ret_Lnk_cells, int64_t *ret_Rnk_cells, int64_t *ret_Tnk_cells, float *ret_Mb)

    # cdef CM_HB_SHADOW_MX *cm_hb_shadow_mx_Create        (CM_t *cm)
    # cdef int              cm_hb_shadow_mx_GrowTo        (CM_t *cm, CM_HB_SHADOW_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit)
    # cdef int              cm_hb_shadow_mx_Dump          (FILE *ofp, CM_t *cm, CM_HB_SHADOW_MX *mx, int print_mx)
    # cdef void             cm_hb_shadow_mx_Destroy       (CM_HB_SHADOW_MX *mx)
    # cdef int              cm_hb_shadow_mx_SizeNeeded    (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int64_t *ret_ny_cells, int64_t *ret_nk_cells, float *ret_Mb)

    # cdef CM_TR_HB_SHADOW_MX *cm_tr_hb_shadow_mx_Create  (CM_t *cm)
    # cdef int              cm_tr_hb_shadow_mx_GrowTo     (CM_t *cm, CM_TR_HB_SHADOW_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit)
    # cdef int              cm_tr_hb_shadow_mx_Dump       (FILE *ofp, CM_t *cm, CM_TR_HB_SHADOW_MX *mx, char mode, int print_mx)
    # cdef void             cm_tr_hb_shadow_mx_Destroy    (CM_TR_HB_SHADOW_MX *mx)
    # cdef int              cm_tr_hb_shadow_mx_SizeNeeded (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int64_t *ret_Jny_cells, int64_t *ret_Lny_cells, int64_t *ret_Rny_cells, 
    #                             int64_t *ret_Jnk_cells, int64_t *ret_Lnk_cells, int64_t *ret_Rnk_cells, int64_t *ret_Tnk_cells, float *ret_Mb)

    # cdef CM_EMIT_MX      *cm_emit_mx_Create     (CM_t *cm)
    # cdef int              cm_emit_mx_GrowTo     (CM_t *cm, CM_EMIT_MX *mx, char *errbuf, int L, float size_limit)
    # cdef int              cm_emit_mx_Dump       (FILE *ofp, CM_t *cm, CM_EMIT_MX *mx, int print_mx)
    # cdef void             cm_emit_mx_Destroy    (CM_EMIT_MX *mx)
    # cdef int              cm_emit_mx_SizeNeeded (CM_t *cm, char *errbuf, int L, int64_t *ret_l_ncells, int64_t *ret_r_ncells, float *ret_Mb)

    # cdef CM_TR_EMIT_MX   *cm_tr_emit_mx_Create     (CM_t *cm)
    # cdef int              cm_tr_emit_mx_GrowTo     (CM_t *cm, CM_TR_EMIT_MX *mx, char *errbuf, int L, float size_limit)
    # cdef int              cm_tr_emit_mx_Dump       (FILE *ofp, CM_t *cm, CM_TR_EMIT_MX *mx, char mode, int print_mx)
    # cdef void             cm_tr_emit_mx_Destroy    (CM_TR_EMIT_MX *mx)
    # cdef int              cm_tr_emit_mx_SizeNeeded (CM_t *cm, char *errbuf, int L, int64_t *ret_l_ncells, int64_t *ret_r_ncells, float *ret_Mb)

    # cdef CM_HB_EMIT_MX   *cm_hb_emit_mx_Create     (CM_t *cm)
    # cdef int              cm_hb_emit_mx_GrowTo     (CM_t *cm, CM_HB_EMIT_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit)
    # cdef int              cm_hb_emit_mx_Dump       (FILE *ofp, CM_t *cm, CM_HB_EMIT_MX *mx, int print_mx)
    # cdef void             cm_hb_emit_mx_Destroy    (CM_HB_EMIT_MX *mx)
    # cdef int              cm_hb_emit_mx_SizeNeeded (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int L, int64_t *ret_l_ncells, int64_t *ret_r_ncells, float *ret_Mb)

    # cdef CM_TR_HB_EMIT_MX *cm_tr_hb_emit_mx_Create     (CM_t *cm)
    # cdef int               cm_tr_hb_emit_mx_GrowTo     (CM_t *cm, CM_TR_HB_EMIT_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit)
    # cdef int               cm_tr_hb_emit_mx_Dump       (FILE *ofp, CM_t *cm, CM_TR_HB_EMIT_MX *mx, char mode, int print_mx)
    # cdef void              cm_tr_hb_emit_mx_Destroy    (CM_TR_HB_EMIT_MX *mx)
    # cdef int               cm_tr_hb_emit_mx_SizeNeeded (CM_t *cm, char *errbf, CP9Bands_t *cp9b, int L, int64_t *ret_l_ncells, int64_t *ret_r_ncells, float *ret_Mb)

    cdef int   cm_scan_mx_Create            (CM_t *cm, char *errbuf, int do_float, int do_int, CM_SCAN_MX **ret_smx)
    cdef int   cm_scan_mx_InitializeFloats  (CM_t *cm, CM_SCAN_MX *smx, char *errbuf)
    cdef int   cm_scan_mx_InitializeIntegers(CM_t *cm, CM_SCAN_MX *smx, char *errbuf)
    cdef float cm_scan_mx_SizeNeeded        (CM_t *cm, int do_float, int do_int)
    cdef void  cm_scan_mx_Destroy           (CM_t *cm, CM_SCAN_MX *smx)
    cdef void  cm_scan_mx_Dump              (FILE *ofp, CM_t *cm, int j, int i0, int qdbidx, int doing_float)

    cdef int   cm_tr_scan_mx_Create            (CM_t *cm, char *errbuf, int do_float, int do_int, CM_TR_SCAN_MX **ret_smx)
    cdef int   cm_tr_scan_mx_InitializeFloats  (CM_t *cm, CM_TR_SCAN_MX *trsmx, char *errbuf)
    cdef int   cm_tr_scan_mx_InitializeIntegers(CM_t *cm, CM_TR_SCAN_MX *trsmx, char *errbuf)
    cdef float cm_tr_scan_mx_SizeNeeded        (CM_t *cm, int do_float, int do_int)
    cdef void  cm_tr_scan_mx_Destroy           (CM_t *cm, CM_TR_SCAN_MX *smx)
    cdef void  cm_tr_scan_mx_Dump              (FILE *ofp, CM_t *cm, int j, int i0, int qdbidx, int doing_float)

    # cdef GammaHitMx_t    *CreateGammaHitMx              (int L, int64_t i0, float cutoff)
    # cdef void             FreeGammaHitMx                (GammaHitMx_t *gamma)
    # cdef int              UpdateGammaHitMx              (CM_t *cm, char *errbuf, int pass_idx, GammaHitMx_t *gamma, int j, int dmin, int dmax, float *bestsc, int *bestr, char *bestmode, int W, double **act)
    # cdef int              ReportHitsGreedily            (CM_t *cm, char *errbuf, int pass_idx, int j, int dmin, int dmax, float *bestsc, int *bestr, char *bestmode, int W, double **act, int64_t i0, int64_tj0, float cutoff, CM_TOPHITS *hitlist)
    # cdef void             TBackGammaHitMx               (GammaHitMx_t *gamma, CM_TOPHITS *hitlist, int64_t i0, int64_t j0)
