cdef extern from "infernal.h" nogil:
    cdef struct expinfo_s:
        double cur_eff_dbsize
        double lambda_ "lambda"
        double mu_extrap
        double mu_orig
        double dbsize
        int    nrandhits
        double tailp
        int    is_valid
    ctypedef expinfo_s ExpInfo_t

    cdef enum:
        EXP_CM_GC
        EXP_CM_GI
        EXP_CM_LC
        EXP_CM_LI
        EXP_NMODES

    # int        debug_print_expinfo_array(CM_t *cm, char *errbuf, ExpInfo_t **expA);
    # int        debug_print_expinfo(ExpInfo_t *exp);
    # int        get_gc_comp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int start, int stop);
    # int        get_alphabet_comp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int start, int stop, float **ret_freq); 
    # int        GetDBSize (ESL_SQFILE *sqfp, char *errbuf, long *ret_N, int *ret_nseq, int *ret_namewidth);
    # int        GetDBInfo(const ESL_ALPHABET *abc, ESL_SQFILE *sqfp, char *errbuf, long *ret_N, int *ret_nseq, double **ret_gc_ct);
    # int        E2ScoreGivenExpInfo(ExpInfo_t *exp, char *errbuf, double E, float *ret_sc);
    # int        P2ScoreGivenExpInfo(ExpInfo_t *exp, char *errbuf, double P, float *ret_sc);
    # double     Score2E(float x, double mu, double lambda, double eff_dbsize);
    # float      cm_p7_E2Score(double E, double Z, int hitlen, float mu, float lambda);
    # float      cm_p7_P2Score(double P, float mu, float lambda);
    # int        ExpModeIsLocal(int exp_mode);
    # int        ExpModeIsInside(int exp_mode);
    # ExpInfo_t *CreateExpInfo();
    # void       SetExpInfo(ExpInfo_t *exp, double lambda, double mu_orig, double dbsize, int nrandhits, double tailp);
    # ExpInfo_t *DuplicateExpInfo(ExpInfo_t *src);
    # char      *DescribeExpMode(int exp_mode);
    # int        UpdateExpsForDBSize(CM_t *cm, char *errbuf, double dbsize);
    # int        CreateGenomicHMM(const ESL_ALPHABET *abc, char *errbuf, double **ret_sA, double ***ret_tAA, double ***ret_eAA, int *ret_nstates);
    # int        SampleGenomicSequenceFromHMM(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, char *errbuf, double *sA, double **tAA, double **eAA, int nstates, int L, ESL_DSQ **ret_dsq);
    # int        CopyExpInfo(ExpInfo_t *src, ExpInfo_t *dest);
