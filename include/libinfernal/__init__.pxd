

cdef extern from "infernal.h" nogil:

    const double DEFAULT_BETA_W
    const double DEFAULT_BETA_QDB1
    const double DEFAULT_BETA_QDB2
    const double DEFAULT_TAU
    const double DEFAULT_PBEGIN
    const double DEFAULT_PEND
    const double DEFAULT_ETARGET
    const double DEFAULT_ETARGET_HMMFILTER
    const double DEFAULT_NULL2_OMEGA
    const double DEFAULT_NULL3_OMEGA
    const double V1P0_NULL2_OMEGA
    const double V1P0_NULL3_OMEGA
    const double DEFAULT_EL_SELFPROB
    const double DEFAULT_MAXTAU
    const double DEFAULT_CP9BANDS_THRESH1
    const double DEFAULT_CP9BANDS_THRESH2
    const double DEFAULT_HB_MXSIZE_MAX_MB
    const double DEFAULT_HB_MXSIZE_MAX_W
    const double DEFAULT_HB_MXSIZE_MIN_MB
    const double DEFAULT_HB_MXSIZE_MIN_W

    const double TAU_MULTIPLIER
    const double MAX_CP9BANDS_THRESH1
    const double MIN_CP9BANDS_THRESH2
    const double DELTA_CP9BANDS_THRESH1
    const double DELTA_CP9BANDS_THRESH2

    const size_t GC_SEGMENTS
    const size_t CM_MAX_RESIDUE_COUNT

    const size_t CM_p7_NEVPARAM
    const float  CM_p7_EVPARAM_UNSET 
    
    cdef enum cm_p7_evparams_e:
        CM_p7_LMMU
        CM_p7_LMLAMBDA
        CM_p7_LVMU
        CM_p7_LVLAMBDA
        CM_p7_LFTAU
        CM_p7_LFLAMBDA
        CM_p7_GFMU
        CM_p7_GFLAMBDA

    const  double IMPOSSIBLE
    const  double MAXSCOREVAL
    const  double IMPROBABLE
    #define NOT_IMPOSSIBLE(x)  ((x) > -9.999e35) 
    #define NOT_IMPROBABLE(x)  ((x) > -4.999e35) 
    #define sreLOG2(x)  ((x) > 0 ? log(x) * 1.44269504 : IMPOSSIBLE)
    #define sreEXP2(x)  (exp((x) * 0.69314718 )) 
    #define epnEXP10(x) (exp((x) * 2.30258509 ))
    #define NOTZERO(x)  (fabs(x - 0.) > -1e6)
    const double INFTY

    const float INTSCALE
    const size_t LOGSUM_TBL

    cdef enum emitmode_e:
        EMITLEFT
        EMITRIGHT
        EMITPAIR
        EMITNONE
    const size_t nEMITMODES


    cdef size_t MAXCONNECT
    cdef enum:
        D_st
        MP_st
        ML_st
        MR_st
        IL_st
        IR_st
        S_st
        E_st
        B_st
        EL_st

    cdef size_t NODETYPES
    cdef enum:
        DUMMY_nd
        BIF_nd
        MATP_nd
        MATL_nd
        MATR_nd
        BEGL_nd		
        BEGR_nd
        ROOT_nd		
        END_nd

    cdef size_t UNIQUESTATES
    cdef enum:
        DUMMY  
        ROOT_S 
        ROOT_IL
        ROOT_IR
        BEGL_S 
        BEGR_S 
        BEGR_IL
        MATP_MP
        MATP_ML
        MATP_MR
        MATP_D
        MATP_IL
        MATP_IR
        MATL_ML
        MATL_D
        MATL_IL
        MATR_MR
        MATR_D
        MATR_IR
        END_E
        BIF_B
        END_EL

    cdef enum:
        TRACE_LEFT_CHILD 
        TRACE_RIGHT_CHILD

    cdef enum:
        PDA_RESIDUE
        PDA_STATE
        PDA_MARKER 