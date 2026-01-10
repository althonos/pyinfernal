from libeasel.alphabet cimport ESL_ALPHABET


cdef extern from "infernal.h" nogil:

    cdef struct cplan9_s:
        const ESL_ALPHABET *abc
        int     M
        float **t
        float **mat
        float **ins

        float  *begin
        float  *end

        int  **tsc
        int  **msc
        int  **isc
        int   *bsc
        int   *esc
        int   *tsc_mem 
        int   *msc_mem
        int   *isc_mem
        int   *bsc_mem
        int   *esc_mem
        int   *otsc

        float  *null
        float  null2_omega
        float  null3_omega
        float  p1

        float  el_self
        int    el_selfsc
        int   *has_el
        int   *el_from_ct
        int  **el_from_idx
        int  **el_from_cmnd
        int flags
    ctypedef cplan9_s CP9_t

    cdef enum:
        CPLAN9_HASBITS
        CPLAN9_HASPROB
        CPLAN9_LOCAL_BEGIN
        CPLAN9_LOCAL_END
        CPLAN9_LOCAL_EL

    cdef enum cp9o_tsc_e:
        cp9O_MM
        cp9O_IM
        cp9O_DM
        cp9O_BM
        cp9O_MI
        cp9O_II
        cp9O_DI
        cp9O_MD
        cp9O_ID
        cp9O_DD
        cp9O_ME
        cp9O_MEL
