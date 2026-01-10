from libc.stdint cimport int64_t, uint32_t, uint64_t
from posix.types cimport off_t

from libeasel cimport eslERRBUFSIZE
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.getopts cimport ESL_GETOPTS
from libeasel.keyhash cimport ESL_KEYHASH
from libeasel.random cimport ESL_RANDOMNESS
from libeasel.sq cimport ESL_SQ
from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_domaindef cimport P7_DOMAINDEF
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_scoredata cimport P7_SCOREDATA
from libinfernal.cm cimport CM_t
from libinfernal.cm_file cimport CM_FILE
from libinfernal.cm_tophits cimport CM_TOPHITS

if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_omx cimport P7_OM_BLOCK
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_omx cimport P7_OM_BLOCK
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon.p7_omx cimport P7_OM_BLOCK
    from libhmmer.impl_neon.p7_oprofile cimport P7_OPROFILE

cdef extern from "infernal.h" nogil:

    cdef enum cm_pipemodes_e:
        CM_SEARCH_SEQS
        CM_SCAN_MODELS

    cdef enum cm_newmodelmodes_e:
        CM_NEWMODEL_MSV
        CM_NEWMODEL_CM

    cdef enum cm_zsetby_e:
        CM_ZSETBY_SSIINFO
        CM_ZSETBY_SSI_AND_QLENGTH
        CM_ZSETBY_FILEREAD
        CM_ZSETBY_OPTION
        CM_ZSETBY_FILEINFO

    cdef enum:
        PLI_PASS_CM_SUMMED
        PLI_PASS_STD_ANY
        PLI_PASS_5P_ONLY_FORCE
        PLI_PASS_3P_ONLY_FORCE
        PLI_PASS_5P_AND_3P_FORCE
        PLI_PASS_5P_AND_3P_ANY
        PLI_PASS_HMM_ONLY_ANY
    const size_t NPLI_PASSES

    cdef struct cm_pipeline_accounting_s:
        uint64_t npli_top
        uint64_t npli_bot
        uint64_t nres_top
        uint64_t nres_bot
        uint64_t n_past_msv
        uint64_t n_past_vit
        uint64_t n_past_fwd
        uint64_t n_past_gfwd
        uint64_t n_past_edef
        uint64_t n_past_cyk
        uint64_t n_past_ins
        uint64_t n_output
        uint64_t n_past_msvbias
        uint64_t n_past_vitbias
        uint64_t n_past_fwdbias
        uint64_t n_past_gfwdbias
        uint64_t n_past_edefbias
        uint64_t pos_past_msv
        uint64_t pos_past_vit
        uint64_t pos_past_fwd
        uint64_t pos_past_gfwd
        uint64_t pos_past_edef
        uint64_t pos_past_cyk
        uint64_t pos_past_ins
        uint64_t pos_output
        uint64_t pos_past_msvbias
        uint64_t pos_past_vitbias
        uint64_t pos_past_fwdbias
        uint64_t pos_past_gfwdbias
        uint64_t pos_past_edefbias
        uint64_t n_overflow_fcyk
        uint64_t n_overflow_final
        uint64_t n_aln_hb
        uint64_t n_aln_dccyk
    ctypedef cm_pipeline_accounting_s CM_PLI_ACCT

    cdef struct cm_pipeline_s:
        # P7_OMX       *oxf
        # P7_OMX       *oxb
        # P7_OMX       *fwd
        # P7_OMX       *bck
        # P7_GMX       *gxf
        # P7_GMX       *gxb
        # P7_GMX       *gfwd
        # P7_GMX       *gbck

        cm_pipemodes_e mode
        ESL_ALPHABET *abc
        CM_FILE      *cmfp
        char[eslERRBUFSIZE] errbuf

        int 		maxW
        int 		cmW
        int 		clen
        int64_t       cur_cm_idx
        int           cur_clan_idx

        int64_t       cur_seq_idx
        int64_t       cur_pass_idx

        uint64_t      nseqs
        uint64_t      nmodels
        uint64_t      nnodes
        uint64_t      nmodels_hmmonly
        uint64_t      nnodes_hmmonly
        CM_PLI_ACCT[NPLI_PASSES] acct

        ESL_RANDOMNESS *r
        bint            do_reseeding
        P7_DOMAINDEF   *ddef

        float         mxsize_limit
        bint          mxsize_set
        bint          be_verbose
        bint          do_top
        bint          do_bot
        bint          show_accessions
        bint          show_alignments
        double        maxtau
        bint          do_wcx
        float         wcx
        bint          do_one_cmpass
        bint          do_one_cmpass_olap
        bint          do_not_iterate
        float         smult
        float         wmult
        float         cmult
        float         mlmult
        bint          do_time_F1
        bint          do_time_F2
        bint          do_time_F3
        bint          do_time_F4
        bint          do_time_F5
        bint          do_time_F6
        bint          do_trm_F3

        bint    by_E
        double  E
        double  T
        int     use_bit_cutoffs

        bint    inc_by_E
        double  incE
        double  incT

        double  Z
        cm_zsetby_e Z_setby

        bint    do_max
        bint    do_nohmm
        bint    do_mid
        bint    do_rfam

        double  F1
        double  F2
        double  F3
        double  F4
        double  F5
        double  F6
        double  F1b
        double  F2b
        double  F3b
        double  F4b
        double  F5b

        bint    do_msv
        bint    do_vit
        bint    do_fwd
        bint    do_gfwd
        bint    do_edef
        bint    do_fcyk
        bint    do_msvbias
        bint    do_vitbias
        bint    do_fwdbias
        bint    do_gfwdbias
        bint    do_edefbias

        bint     do_trunc_ends
        bint     do_trunc_any
        bint     do_trunc_int
        bint     do_trunc_only
        bint     do_trunc_5p_ends
        bint     do_trunc_3p_ends

        float  rt1
        float  rt2
        float  rt3
        int    ns

        bint    do_fcykenv
        double  F6env
        bint    do_null2
        bint    do_null3
        bint    do_glocal_cm_always
        bint    do_glocal_cm_cur
        bint    do_glocal_cm_sometimes
        int     fcyk_cm_search_opts
        int     final_cm_search_opts
        int     fcyk_cm_exp_mode
        int     final_cm_exp_mode
        double  fcyk_beta
        double  final_beta
        double  fcyk_tau
        double  final_tau

        bint    do_hmmonly_cur
        bint    do_hmmonly_always
        bint    do_hmmonly_never
        bint    do_max_hmmonly

        double  F1_hmmonly
        double  F2_hmmonly
        double  F3_hmmonly

        bint    do_bias_hmmonly
        bint    do_null2_hmmonly

        int     cm_config_opts
        int     cm_align_opts
    ctypedef cm_pipeline_s CM_PIPELINE

    CM_PIPELINE *cm_pipeline_Create (ESL_GETOPTS *go, ESL_ALPHABET *abc, int clen_hint, int L_hint, int64_t Z, cm_zsetby_e Z_setby, cm_pipemodes_e mode)
    int          cm_pipeline_Reuse  (CM_PIPELINE *pli)
    void         cm_pipeline_Destroy(CM_PIPELINE *pli, CM_t *cm)
    int          cm_pipeline_Merge  (CM_PIPELINE *p1, CM_PIPELINE *p2)

    int   cm_pli_TargetReportable  (CM_PIPELINE *pli, float score,     double Eval)
    int   cm_pli_TargetIncludable  (CM_PIPELINE *pli, float score,     double Eval)
    int   cm_pli_NewModel          (CM_PIPELINE *pli, int modmode, CM_t *cm, int cm_clen, int cm_W, int cm_nbp, P7_OPROFILE *om, P7_BG *bg, float *p7_evparam, int p7_max_length, int64_t cur_cm_idx, int cur_clan_idx, ESL_KEYHASH *glocal_kh)
    int   cm_pli_NewModelThresholds(CM_PIPELINE *pli, CM_t *cm)
    int   cm_pli_NewSeq            (CM_PIPELINE *pli, const ESL_SQ *sq, int64_t cur_seq_idx)
    int   cm_Pipeline              (CM_PIPELINE *pli, off_t cm_offset, P7_OPROFILE *om, P7_BG *bg, float *p7_evparam, P7_SCOREDATA *msvdata, ESL_SQ *sq, CM_TOPHITS *hitlist, int in_rc, P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_PROFILE **opt_Rgm, P7_PROFILE **opt_Lgm, P7_PROFILE **opt_Tgm, CM_t **opt_cm)
    # int   cm_pli_Statistics    (FILE *ofp, CM_PIPELINE *pli, ESL_STOPWATCH *w)
    int   cm_pli_ZeroAccounting(CM_PLI_ACCT *pli_acct)
    int   cm_pli_PassEnforcesFirstRes(int pass_idx)
    int   cm_pli_PassEnforcesFinalRes(int pass_idx)
    int   cm_pli_PassAllowsTruncation(int pass_idx)
    void  cm_pli_AdjustNresForOverlaps(CM_PIPELINE *pli, int64_t noverlap, int in_rc)
