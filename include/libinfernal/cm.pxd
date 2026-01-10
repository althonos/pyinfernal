from libc.stdint cimport uint32_t, uint64_t
from posix.types cimport off_t

from libeasel.alphabet cimport ESL_ALPHABET
from libhmmer.p7_hmm cimport P7_HMM
from libinfernal cimport CM_p7_NEVPARAM
from libinfernal.cm_qdband cimport CM_QDBINFO
from libinfernal.stats cimport ExpInfo_t


cdef extern from "infernal.h" nogil:

    cdef enum:
        CMH_BITS
        CMH_ACC
        CMH_DESC
        CMH_RF
        CMH_GA
        CMH_TC
        CMH_NC
        CMH_CHKSUM
        CMH_MAP
        CMH_CONS
        CMH_LOCAL_BEGIN
        CMH_LOCAL_END
        CMH_EXPTAIL_STATS
        CMH_CP9
        CMH_CP9_TRUNC
        CMH_MLP7
        CMH_FP7
        CM_IS_SUB
        CM_IS_RSEARCH
        CM_RSEARCHTRANS
        CM_RSEARCHEMIT
        CM_EMIT_NO_LOCAL_BEGINS
        CM_EMIT_NO_LOCAL_ENDS
        CM_IS_CONFIGURED

    cdef enum:
        CM_CONFIG_LOCAL
        CM_CONFIG_HMMLOCAL
        CM_CONFIG_HMMEL
        CM_CONFIG_QDB
        CM_CONFIG_W_BETA
        CM_CONFIG_TRUNC
        CM_CONFIG_SCANMX
        CM_CONFIG_TRSCANMX
        CM_CONFIG_SUB
        CM_CONFIG_NONBANDEDMX

    cdef enum:
        CM_ALIGN_HBANDED
        CM_ALIGN_P7BANDED
        CM_ALIGN_NONBANDED
        CM_ALIGN_CYK
        CM_ALIGN_OPTACC
        CM_ALIGN_SAMPLE
        CM_ALIGN_POST
        CM_ALIGN_SMALL
        CM_ALIGN_SUMS
        CM_ALIGN_SUB
        CM_ALIGN_HMMVITERBI
        CM_ALIGN_CHECKINOUT
        CM_ALIGN_CHECKPARSESC
        CM_ALIGN_PRINTTREES
        CM_ALIGN_HMMSAFE
        CM_ALIGN_SCOREONLY
        CM_ALIGN_FLUSHINSERTS
        CM_ALIGN_CHECKFB
        CM_ALIGN_HMM2IJOLD
        CM_ALIGN_QDB
        CM_ALIGN_INSIDE
        CM_ALIGN_TRUNC
        CM_ALIGN_XTAU

    cdef enum:
        CM_SEARCH_HBANDED
        CM_SEARCH_QDB
        CM_SEARCH_NONBANDED
        CM_SEARCH_HMMALNBANDS
        CM_SEARCH_SUMS
        CM_SEARCH_INSIDE
        CM_SEARCH_NOALIGN
        CM_SEARCH_RSEARCH
        CM_SEARCH_CMNOTGREEDY
        CM_SEARCH_HMM2IJOLD
        CM_SEARCH_NULL2
        CM_SEARCH_NULL3

    cdef enum cm_w_setby_e:
        CM_W_SETBY_INIT
        CM_W_SETBY_CMFILE
        CM_W_SETBY_BANDCALC
        CM_W_SETBY_CMDLINE
        CM_W_SETBY_SUBCOPY

    cdef struct cm_s:
        char    *name
        char    *acc
        char    *desc
        char    *rf
        char    *consensus
        uint32_t checksum
        int     *map

        # char  *comlog
        # char  *ctime
        # int    nseq
        # float  eff_nseq
        # float  ga
        # float  tc
        # float  nc

        # float *null

        int   M
        int   clen
        char *sttype
        int  *ndidx
        char *stid

        # int  *cfirst
        # int  *cnum
        # int  *plast
        # int  *pnum

        # int   nodes
        # int  *nodemap
        # char *ndtype

        # float **t
        # float **e
        # float  *begin
        # float  *end

        # float **tsc
        # float **esc
        # float **oesc
        # float **lmesc
        # float **rmesc
        # float *beginsc
        # float *endsc

        # int  **itsc
        # int  **iesc
        # int  **ioesc
        # int  **ilmesc
        # int  **irmesc
        # int   *ibeginsc
        # int   *iendsc

        # float  pbegin
        # float  pend

        # float  null2_omega
        # float  null3_omega

        int    flags

        int    W
        double beta_W
        cm_w_setby_e W_setby

        CM_QDBINFO *qdbinfo

        double  tau
        double  maxtau

        int         config_opts
        int         align_opts
        int         search_opts
        float*      root_trans

        float  el_selfsc
        int   iel_selfsc

        # CP9_t      *cp9
        # CP9_t      *Lcp9
        # CP9_t      *Rcp9
        # CP9_t      *Tcp9
        # CP9Map_t   *cp9map
        # CP9Bands_t *cp9b

        # CM_HB_MX           *hb_mx
        # CM_HB_MX           *hb_omx
        # CM_HB_EMIT_MX      *hb_emx
        # CM_HB_SHADOW_MX    *hb_shmx

        # CM_TR_HB_MX        *trhb_mx
        # CM_TR_HB_MX        *trhb_omx
        # CM_TR_HB_EMIT_MX   *trhb_emx
        # CM_TR_HB_SHADOW_MX *trhb_shmx

        # CM_MX              *nb_mx
        # CM_MX              *nb_omx
        # CM_EMIT_MX         *nb_emx
        # CM_SHADOW_MX       *nb_shmx

        # CM_TR_MX           *trnb_mx
        # CM_TR_MX           *trnb_omx
        # CM_TR_EMIT_MX      *trnb_emx
        # CM_TR_SHADOW_MX    *trnb_shmx

        # CM_SCAN_MX         *smx
        # CM_TR_SCAN_MX      *trsmx
        # CP9_MX             *cp9_mx
        # CP9_MX             *cp9_bmx

        ExpInfo_t       **expA

        P7_HMM *mlp7
        P7_HMM *fp7
        float[CM_p7_NEVPARAM] fp7_evparam

        const  ESL_ALPHABET *abc
        off_t  offset

        # CMEmitMap_t     *emap
        # CMConsensus_t   *cmcons
        # CM_TR_PENALTIES *trp
    ctypedef cm_s CM_t

    CM_t *CreateCM(int nnodes, int nstates, int clen, const ESL_ALPHABET *abc)
    CM_t *CreateCMShell()
    void  CreateCMBody(CM_t *cm, int nnodes, int nstates, int clen, const ESL_ALPHABET *abc)
    void  CMZero(CM_t *cm)
    CMRenormalize(CM_t *cm)
    void  FreeCM(CM_t *cm)
    void  CMSimpleProbify(CM_t *cm)
    # int   rsearch_CMProbifyEmissions(CM_t *cm, fullmat_t *fullmat);
    int   CMLogoddsify(CM_t *cm)
    int   CMCountStatetype(CM_t *cm, char type)
    int   CMCountNodetype(CM_t *cm, char type)
    # int   CMSegmentCountStatetype(CM_t *cm, int r, int z, char type);
    # int   CMSubtreeCountStatetype(CM_t *cm, int v, char type);
    # int   CMSubtreeCountNodetype(CM_t *cm, int v, char type);
    # int   CMSubtreeFindEnd(CM_t *cm, int v);
    # int   CalculateStateIndex(CM_t *cm, int node, char utype);
    # int   TotalStatesInNode(int ndtype);
    # int   SplitStatesInNode(int ndtype);
    # int   InsertStatesInNode(int ndtype);
    # int   StateDelta(int sttype);
    # int   StateLeftDelta(int sttype);
    # int   StateRightDelta(int sttype);
    # int   Emitmode(int sttype);
    # int   NumReachableInserts(int stid);
    # void  PrintCM(FILE *fp, CM_t *cm);
    # void  SummarizeCM(FILE *fp, CM_t *cm);
    char *Statetype(int type)
    int   StateCode(char *s)
    char *Nodetype(int type)
    int   NodeCode(char *s)
    char *UniqueStatetype(int type)
    int   UniqueStateCode(char *s)
    int   DeriveUniqueStateCode(int ndtype, int sttype)
    int   StateMapsLeft(char st)
    int   StateMapsRight(char st)
    int   StateMapsMatch(char st)
    int   StateMapsInsert(char st)
    int   StateMapsDelete(char st)
    int   NodeMapsLeft(char ndtype)
    int   NodeMapsRight(char ndtype)
    int   StateIsDetached(CM_t *cm, int v)
    int   CMRebalance(CM_t *cm, char *errbuf, CM_t **ret_new_cm)
    int **IMX2Alloc(int rows, int cols)
    void  IMX2Free(int **mx)
    # float rsearch_calculate_gap_penalty (char from_state, char to_state, int from_node, int to_node, float input_alpha, float input_beta, float input_alphap, float input_betap);
    # int   cm_Exponentiate(CM_t *cm, double z);
    # int   cm_p7_Exponentiate(P7_HMM *hmm, double z);
    # void  cm_banner(FILE *fp, char *progname, char *banner);
    # void  cm_CalcExpSc(CM_t *cm, float **ret_expsc, float **ret_expsc_noss);
    # int   cm_Validate(CM_t *cm, float tol, char *errbuf);
    # char *CMStatetype(char st);
    # char *CMNodetype(char nd);
    # char *CMStateid(char st);
    # char *MarginalMode(char mode);
    # int   ModeEmitsLeft(char mode);
    # int   ModeEmitsRight(char mode);
    # int   cm_SetName(CM_t *cm, char *name);
    # int   cm_SetAccession(CM_t *cm, char *acc);
    # int   cm_SetDescription(CM_t *cm, char *desc);
    # int   cm_SetConsensus(CM_t *cm, CMConsensus_t *cons, ESL_SQ *sq);
    # int   cm_AppendComlog(CM_t *cm, int argc, char **argv, int add_seed, uint32_t seed);
    # int   cm_SetCtime(CM_t *cm);
    # int   DefaultNullModel(const ESL_ALPHABET *abc, float **ret_null);
    # int   CMAllocNullModel(CM_t *cm);
    # void  CMSetNullModel(CM_t *cm, float *null);
    # int   CMReadNullModel(const ESL_ALPHABET *abc, char *nullfile, float **ret_null);
    # int   IntMaxDigits();
    # int   IntDigits(int i);
    # int        cm_GetAvgHitLen(CM_t *cm, char *errbuf, float *ret_avgL_loc, float *ret_avgL_glb);
    # int        CompareCMGuideTrees(CM_t *cm1, CM_t *cm2);
    # void       DumpCMFlags(FILE *fp, CM_t *cm);
    # ESL_GETOPTS *cm_CreateDefaultApp(ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage);
    # CM_P7_OM_BLOCK *cm_p7_oprofile_CreateBlock(int size);
    # void            cm_p7_oprofile_DestroyBlock(CM_P7_OM_BLOCK *block);
    # float **FCalcOptimizedEmitScores      (CM_t *cm);
    # int   **ICalcOptimizedEmitScores      (CM_t *cm);
    # int   **ICopyOptimizedEmitScoresFromFloats(CM_t *cm, float **oesc);
    # int     CloneOptimizedEmitScores      (const CM_t *src, CM_t *dest, char *errbuf);
    # void    DumpOptimizedEmitScores       (CM_t *cm, FILE *fp);
    # void    FreeOptimizedEmitScores       (float **fesc_vAA, int **iesc_vAA, int M);
    # float **FCalcInitDPScores             (CM_t *cm);
    # int   **ICalcInitDPScores             (CM_t *cm);
    # int     cm_nonconfigured_Verify(CM_t *cm, char *errbuf);
    int     cm_Clone(CM_t *cm, char *errbuf, CM_t **ret_cm)
    float   cm_Sizeof(CM_t *cm)
    int     Prob2Score(float p, float null)
    float   Score2Prob(int sc, float null)
    float   Scorify(int sc)
    double *cm_ExpectedStateOccupancy(CM_t *cm)
    # int     cm_ExpectedPositionOccupancy(CM_t *cm, float **ret_mexpocc, float **ret_iexpocc, double **opt_psi, int **opt_m2v_1, int **opt_m2v_2, int **opt_i2v);
    char ***cm_CreateTransitionMap()
    void    cm_FreeTransitionMap(char ***tmap)
    # void    InsertsGivenNodeIndex(CM_t *cm, int nd, int *ret_i1, int *ret_2);
    # int     cm_Guidetree(CM_t *cm, char *errbuf, ESL_MSA *msa, Parsetree_t **ret_gtr);