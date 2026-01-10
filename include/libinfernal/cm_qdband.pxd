cdef extern from "infernal.h" nogil:

    cdef struct cm_s:
        pass
    ctypedef cm_s CM_t

    cdef enum cm_qdbinfo_setby_e:
        CM_QDBINFO_SETBY_INIT
        CM_QDBINFO_SETBY_CMFILE
        CM_QDBINFO_SETBY_BANDCALC
        CM_QDBINFO_SETBY_SUBINIT

    cdef struct cm_qdbinfo_s:
        int M
        double beta1
        int* dmin1
        int* dmax1
        double beta2
        int* dmin2
        int* dmax2
        cm_qdbinfo_setby_e setby
    ctypedef cm_qdbinfo_s CM_QDBINFO

    void     BandExperiment(CM_t *cm)
    # int      CalculateQueryDependentBands(CM_t *cm, char *errbuf, CM_QDBINFO *qdbinfo, double beta_W, int *ret_W, double **ret_gamma0_loc, double **ret_gamma0_glb, int *ret_Z);
    # int      BandCalculationEngine(CM_t *cm, int Z, CM_QDBINFO *qdbinfo, double beta_W, int *ret_W, double ***ret_gamma, double **ret_gamma0_loc, double **ret_gamma0_glb);
    # int      BandTruncationNegligible(double *density, int b, int Z, double *ret_beta);
    # int      BandMonteCarlo(CM_t *cm, int nsample, int Z, double ***ret_gamma);
    # void     FreeBandDensities(CM_t *cm, double **gamma);
    # void     BandBounds(double **gamma, int M, int Z, double p, 
    #             int **ret_min, int **ret_max);
    # void     PrintBandGraph(FILE *fp, double **gamma, int *min, int *max, int v, int Z);
    # void     PrintDPCellsSaved(CM_t *cm, int *min, int *max, int W);
    # void     ExpandBands(CM_t *cm, int qlen, int *dmin, int *dmax);
    # void     qdb_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *dmin, int *dmax, int bdump_level);
    # CM_QDBINFO  *CreateCMQDBInfo(int M, int clen);
    # float        SizeofCMQDBInfo(CM_QDBINFO *qdbinfo);
    # void         FreeCMQDBInfo(CM_QDBINFO *qdbinfo);
    # int          CopyCMQDBInfo(const CM_QDBINFO *src, CM_QDBINFO *dst, char *errbuf);
    # void         DumpCMQDBInfo(FILE *fp, CM_t *cm, CM_QDBINFO *qdbinfo);
    int          CheckCMQDBInfo(CM_QDBINFO *qdbinfo, double beta1, bint do_check1, double beta2, bint do_check2)
