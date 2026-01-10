from libc.stdio cimport FILE
from posix.types cimport off_t

from libeasel cimport eslERRBUFSIZE
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.fileparser cimport ESL_FILEPARSER
from libeasel.ssi cimport ESL_SSI
from libhmmer.p7_hmmfile cimport P7_HMMFILE
from libinfernal.cm cimport CM_t


cdef extern from "infernal.h" nogil:

    cdef enum cm_file_formats_e:
        CM_FILE_1
        CM_FILE_1a

    cdef struct cm_file_s:
        FILE         *f
        char         *fname
        ESL_SSI      *ssi

        bint           do_gzip
        bint           do_stdin
        bint           newly_opened
        bint           is_binary
        bint           is_pressed

        int            format
        int           (*parser)(cm_file_s*, int, ESL_ALPHABET**, CM_t**) 
        ESL_FILEPARSER *efp

        P7_HMMFILE   *hfp

        FILE         *ffp
        FILE         *pfp
        char          errbuf[eslERRBUFSIZE]
        char          msv_errbuf[eslERRBUFSIZE]
    ctypedef cm_file_s CM_FILE

    int read_asc_1p1_cm(CM_FILE *hfp, int read_fp7, ESL_ALPHABET **ret_abc, CM_t **opt_cm)
    int read_bin_1p1_cm(CM_FILE *hfp, int read_fp7, ESL_ALPHABET **ret_abc, CM_t **opt_cm)
    int read_asc_1p0_cm(CM_FILE *hfp, int read_fp7, ESL_ALPHABET **ret_abc, CM_t **opt_cm)

    int     cm_file_Open(char *filename, char *env, bint allow_1p0, CM_FILE **ret_cmfp, char *errbuf) except *
    int     cm_file_OpenNoDB(char *filename, char *env, bint allow_1p0, CM_FILE **ret_cmfp, char *errbuf) except *
    int     cm_file_OpenBuffer(char *buffer, int size, bint allow_1p0, CM_FILE **ret_cmfp) except *
    void    cm_file_Close(CM_FILE *cmfp)
    int     cm_file_CreateLock(CM_FILE *cmfp)
    int     cm_file_WriteASCII(FILE *fp, int format, CM_t *cm) except *
    int     cm_file_WriteBinary(FILE *fp, int format, CM_t *cm, off_t *opt_fp7_offset) except *
    int     cm_file_Read(CM_FILE *cmfp, int read_fp7, ESL_ALPHABET **ret_abc, CM_t **opt_cm) except *
    int     cm_file_PositionByKey(CM_FILE *cmfp, const char *key)
    int     cm_file_Position(CM_FILE *cmfp, const off_t offset)
    # int     cm_p7_hmmfile_Read(CM_FILE *cmfp, ESL_ALPHABET *abc, off_t offset, P7_HMM **ret_hmm)
    # int     cm_p7_oprofile_Write(FILE *ffp, FILE *pfp, off_t cm_offset, int cm_len, int cm_W, int cm_nbp, float gfmu, float gflambda, P7_OPROFILE *om)
    # int     cm_p7_oprofile_ReadMSV(CM_FILE *cmfp, int read_scores, ESL_ALPHABET **byp_abc, off_t *ret_cm_offset, int *ret_cm_clen, int *ret_cm_W, int *ret_cm_nbp, float *ret_gfmu, float *ret_gflambda, P7_OPROFILE **ret_om)
    # int     cm_p7_oprofile_ReadBlockMSV(CM_FILE *cmfp, int64_t cm_idx, ESL_ALPHABET **byp_abc, CM_P7_OM_BLOCK *hmmBlock)
    # int     cm_p7_oprofile_Position(CM_FILE *cmfp, off_t offset)
    # int     cm_file_Write1p0ASCII(FILE *fp, CM_t *cm)
