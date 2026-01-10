from libc cimport errno
from libc.stdio cimport FILE, fopen, fclose
from libc.stdint cimport uint32_t, uint64_t, int64_t
from libc.stdlib cimport malloc, calloc, free
from libc.string cimport memset, memcpy, memmove, strdup, strlen

cimport libeasel
cimport libeasel.vec
cimport libeasel.fileparser
cimport libhmmer.p7_bg
cimport libhmmer.p7_profile
cimport libhmmer.p7_hmm
cimport libhmmer.p7_scoredata
cimport libhmmer.modelconfig
cimport libinfernal.cm
cimport libinfernal.cm_mx
cimport libinfernal.cm_file
cimport libinfernal.cm_tophits
cimport libinfernal.cm_pipeline
cimport libinfernal.cm_qdband
cimport libinfernal.cm_modelconfig
cimport libinfernal.cm_p7_modelconfig
from libeasel cimport eslERRBUFSIZE
from libeasel.alphabet cimport ESL_ALPHABET
from libeasel.fileparser cimport ESL_FILEPARSER
from libeasel.sq cimport ESL_SQ
from libeasel.random cimport ESL_RANDOMNESS
from libhmmer.p7_hmm cimport P7_HMM
from libhmmer.p7_hmmfile cimport P7_HMMFILE
from libhmmer.logsum cimport p7_FLogsumInit
from libinfernal cimport CM_p7_NEVPARAM
from libinfernal.cm_file cimport CM_FILE, cm_file_formats_e
from libinfernal.cm_pipeline cimport CM_PIPELINE, cm_zsetby_e, cm_pipemodes_e, cm_newmodelmodes_e
from libinfernal.cm_tophits cimport CM_TOPHITS, CM_HIT
from libinfernal.cm cimport CM_t
from libinfernal.cmsearch cimport WORKER_INFO
from libinfernal.logsum cimport FLogsumInit, init_ilogsum

if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_omx cimport P7_OM_BLOCK
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE, p7_oprofile_Create, p7_oprofile_Convert
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_omx cimport P7_OM_BLOCK
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE, p7_oprofile_Create, p7_oprofile_Convert
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon.p7_omx cimport P7_OM_BLOCK
    from libhmmer.impl_neon.p7_oprofile cimport P7_OPROFILE, p7_oprofile_Create, p7_oprofile_Convert

if TARGET_SYSTEM == "Linux":
    from .fileobj.linux cimport fileobj_linux_open as fopen_obj
elif TARGET_SYSTEM == "Darwin" or TARGET_SYSTEM.endswith("BSD"):
    from .fileobj.bsd cimport fileobj_bsd_open as fopen_obj

from pyhmmer.easel cimport (
    Alphabet,
    DigitalSequenceBlock,
    Randomness,
    SequenceFile,
)
from pyhmmer.plan7 cimport (
    HMM,
    Profile,
    OptimizedProfile,
)

include "exceptions.pxi"
# include "_getid.pxi"

# --- Python imports ---------------------------------------------------------

import io
import os
import sys
import warnings
from pyhmmer.errors import (
    UnexpectedError,
    AllocationError,
    AlphabetMismatch,
    EaselError,
    InvalidParameter
)

# --- Constants --------------------------------------------------------------

__version__ = PROJECT_VERSION

cdef dict CM_FILE_FORMATS = {
    "2.0": cm_file_formats_e.CM_FILE_1,
    "3/a": cm_file_formats_e.CM_FILE_1a,
}

cdef dict CM_FILE_MAGIC = {
    # v1a_magic:  cm_file_formats_e.CM_FILE_1,
    # v1a_fmagic: cm_file_formats_e.CM_FILE_1a,
    0xe3edb0b2: cm_file_formats_e.CM_FILE_1,
    0xb1e1e6f3: cm_file_formats_e.CM_FILE_1a,
}

# --- Fused types ------------------------------------------------------------

ctypedef fused SearchTargets:
    SequenceFile
    DigitalSequenceBlock

# --- Cython classes ---------------------------------------------------------

cdef class CM:
    """A data structure storing an Infernal Covariance Model.
    """
    cdef CM_t*              _cm
    cdef readonly Alphabet alphabet
    cdef readonly HMM      filter_hmm
    cdef readonly HMM      ml_hmm

    @staticmethod
    cdef from_ptr(CM_t* cm, Alphabet alphabet = None):
        cdef CM obj = CM.__new__(CM)
        obj._cm = cm

        if alphabet is None:
            obj.alphabet = Alphabet.__new__(Alphabet)
            obj.alphabet._abc = cm.abc
        else:
            obj.alphabet = alphabet

        if cm.fp7 is not NULL:
            obj.filter_hmm = HMM.__new__(HMM)
            obj.filter_hmm._hmm = cm.fp7
            obj.filter_hmm.alphabet = obj.alphabet

        if cm.mlp7 is not NULL:
            obj.ml_hmm = HMM.__new__(HMM)
            obj.ml_hmm._hmm = cm.mlp7
            obj.ml_hmm.alphabet = obj.alphabet

        return obj

    def __cinit__(self):
        self.alphabet = None
        self.filter_hmm = None
        self.ml_hmm = None
        self._cm = NULL

    def __dealloc__(self):
        if self._cm is not NULL:
            self._cm.fp7 = NULL # owned by `self.filter_hmm`
            self._cm.mlp7 = NULL # owned by `self.ml_hmm`
            libinfernal.cm.FreeCM(self._cm)

    @property
    def M(self):
        assert self._cm != NULL
        return self._cm.M

    @property
    def name(self):
        assert self._cm != NULL
        return <bytes> self._cm.name

    @property
    def accession(self):
        """`bytes` or `None`: The accession of the CM, if any.
        """
        assert self._cm != NULL
        return None if self._cm.acc == NULL else <bytes> self._cm.acc


cdef class CMFile:
    """A wrapper around a file storing serialized CMs.

    Example:
        Load the first CM from a CM file located on the
        local filesystem::

            >>> with CMFile("tests/data/cms/iss33.cm") as cm_file:
            ...     cm = cm_file.read()
            >>> cm.name
            b'LSU_rRNA_bacteria'
            >>> cm.accession
            b'RF02541'

        Load all the CMs from a CM file into a `list`::

            >>> with CMFile("tests/data/cms/5.c.cm") as cm_file:
            ...     cms = list(cm_file)
            >>> len(cms)
            5
            >>> [cm.accession for cm in cms]
            [b'RF00005', b'RF00006', b'RF01185', b'RF01855', b'RF01852']

    """

    cdef CM_FILE* _fp
    cdef str      _name
    cdef Alphabet _alphabet

    # --- Constructor --------------------------------------------------------

    @staticmethod
    cdef CM_FILE* _open_fileobj(object fh) except *:
        cdef int         status
        cdef char*       token
        cdef int         token_len
        cdef bytes       filename
        cdef object      fh_       = fh
        cdef CM_FILE*    cmfp      = NULL

        # use buffered IO to be able to peek efficiently
        if not hasattr(fh, "peek"):
            fh_ = io.BufferedReader(fh)

        # attempt to allocate space for the P7_HMMFILE
        cmfp = <CM_FILE*> malloc(sizeof(CM_FILE))
        if cmfp == NULL:
            raise AllocationError("CM_FILE", sizeof(CM_FILE))

        # store options
        cmfp.f            = fopen_obj(fh_, "r")
        cmfp.do_stdin     = False
        cmfp.do_gzip      = True
        cmfp.newly_opened = True
        cmfp.is_pressed   = False
        cmfp.is_binary    = False

        # set pointers as NULL for now
        cmfp.parser    = NULL
        cmfp.efp       = NULL
        cmfp.ffp       = NULL
        cmfp.hfp       = NULL
        cmfp.pfp       = NULL
        cmfp.ssi       = NULL
        cmfp.fname     = NULL
        cmfp.errbuf[0] = b"\0"

        # set up the HMM file 
        cmfp.hfp = <P7_HMMFILE*> malloc(sizeof(P7_HMMFILE))
        if cmfp.hfp == NULL:
            libinfernal.cm_file.cm_file_Close(cmfp)
            raise AllocationError("P7_HMMFILE", sizeof(P7_HMMFILE))
        cmfp.hfp.do_gzip      = cmfp.do_gzip
        cmfp.hfp.do_stdin     = cmfp.do_stdin
        cmfp.hfp.newly_opened = True
        cmfp.hfp.is_pressed   = cmfp.is_pressed
        cmfp.hfp.parser       = NULL
        cmfp.hfp.efp          = NULL
        cmfp.hfp.ffp          = NULL
        cmfp.hfp.pfp          = NULL
        cmfp.hfp.ssi          = NULL
        cmfp.hfp.fname        = NULL
        cmfp.hfp.errbuf[0]    = '\0'

        # NOTE(@althonos): Because we set `do_gzip=True`, the parser will now
        #                  expect a lot of things to be available only through
        #                  streams, and won't attempt to e.g. `seek` the 
        #                  internal file object (or at least not as often).
        cmfp.hfp.f            = cmfp.f

        # extract the filename if the file handle has a `name` attribute
        if getattr(fh, "name", None) is not None:
            filename = fh.name.encode()
            cmfp.fname = strdup(filename)
            if cmfp.fname == NULL:
                libinfernal.cm_file.cm_file_Close(cmfp)
                raise AllocationError("char", sizeof(char), strlen(filename))
            cmfp.hfp.fname = strdup(filename)
            if cmfp.hfp.fname == NULL:
                libinfernal.cm_file.cm_file_Close(cmfp)
                raise AllocationError("char", sizeof(char), strlen(filename))

        # check if the parser is in binary format,
        magic = int.from_bytes(fh_.peek(4)[:4], sys.byteorder)
        if magic in CM_FILE_MAGIC:
            cmfp.format = CM_FILE_MAGIC[magic]
            cmfp.parser = libinfernal.cm_file.read_bin_1p1_cm
            cmfp.is_binary = True
            # NB: the file must be advanced, since read_bin_1p1_cm assumes
            #     the binary tag has been skipped already, buf we only peeked
            #     so far; note that we advance without seeking or rewinding.
            fh_.read(4)
            return cmfp
        elif (magic & 0x80000000) != 0:
            raise ValueError(f"Format tag appears binary, but unrecognized: 0x{magic:08x}")

        # create and configure the file parser
        cmfp.efp = libeasel.fileparser.esl_fileparser_Create(cmfp.f)
        if cmfp.efp == NULL:
            libinfernal.cm_file.cm_file_Close(cmfp)
            raise AllocationError("ESL_FILEPARSER", sizeof(ESL_FILEPARSER))
        status = libeasel.fileparser.esl_fileparser_SetCommentChar(cmfp.efp, b"#")
        if status != libeasel.eslOK:
            libinfernal.cm_file.cm_file_Close(cmfp)
            raise UnexpectedError(status, "esl_fileparser_SetCommentChar")

        # get the magic string at the beginning
        status = libeasel.fileparser.esl_fileparser_NextLine(cmfp.efp)
        if status == libeasel.eslEOF:
            raise EOFError("CM file is empty")
        elif status != libeasel.eslOK:
            libinfernal.cm_file.cm_file_Close(cmfp)
            raise UnexpectedError(status, "esl_fileparser_NextLine")
        status = libeasel.fileparser.esl_fileparser_GetToken(cmfp.efp, &token, &token_len)
        if status != libeasel.eslOK:
            libinfernal.cm_file.cm_file_Close(cmfp)
            raise UnexpectedError(status, "esl_fileparser_GetToken")

        # detect the format
        if token == b"INFERNAL1/a":
            cmfp.parser = libinfernal.cm_file.read_asc_1p1_cm
            cmfp.format = cm_file_formats_e.CM_FILE_1a
        elif token == b"INFERNAL-1":
            cmfp.parser = libinfernal.cm_file.read_asc_1p0_cm
            cmfp.format = cm_file_formats_e.CM_FILE_1

        # check the format tag was recognized
        if cmfp.parser == NULL:
            text = token.decode("utf-8", "replace")
            libinfernal.cm_file.cm_file_Close(cmfp)
            raise ValueError("Unrecognized format tag in CM file: {!r}".format(text))

        # return the finalized CM_FILE*
        return cmfp

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._alphabet = None
        self._fp = NULL
        self._name = None

    def __init__(self, object file, bint db = True):
        cdef int                 status
        cdef bytes               fspath
        cdef char[eslERRBUFSIZE] errbuf

        try:
            fspath = os.fsencode(file)
            self._name = os.fsdecode(fspath)
            if db:
                function = "cm_file_Open"
                status = libinfernal.cm_file.cm_file_Open(fspath, NULL, True, &self._fp, errbuf)
            else:
                function = "cm_file_OpenNoDb"
                status = libinfernal.cm_file.cm_file_OpenNoDB(fspath, NULL, True, &self._fp, errbuf)
        except TypeError as e:
            self._fp = CMFile._open_fileobj(file)
            status   = libeasel.eslOK

        if status == libeasel.eslENOTFOUND:
            raise FileNotFoundError(errno.ENOENT, f"No such file or directory: {file!r}")
        elif status == libeasel.eslEFORMAT:
            if fspath is not None:
                if os.path.isdir(fspath):
                    raise IsADirectoryError(errno.EISDIR, f"Is a directory: {file!r}")
                elif os.stat(file).st_size == 0:
                    raise EOFError("CM file is empty")
            raise ValueError("format not recognized by Infernal")
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, function)

        self._alphabet = Alphabet.__new__(Alphabet)
        self._alphabet._abc = NULL

    def __dealloc__(self):
        if self._fp:
            warnings.warn("unclosed CM file", ResourceWarning)
            self.close()

    def __repr__(self):
        cdef str ty = type(self).__name__
        if self._name is not None:
            return f"{ty}({self._name!r})"
        else:
            return super().__repr__()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        cdef CM cm = self.read()
        if cm is None:
            raise StopIteration()
        return cm

    # --- Properties ---------------------------------------------------------

    @property
    def closed(self):
        """`bool`: Whether the `CMFile` is closed or not.
        """
        return self._fp == NULL

    @property
    def name(self):
        """`str` or `None`: The path to the CM file, if known.
        """
        return self._name

    # --- Python Methods -----------------------------------------------------

    cpdef CM read(self):
        """Read the next CM from the file.

        Returns:
            `~pyinfernal.cm.CM` or `None`: The next CM in the file, or
            `None` if all CMs were read from the file already.

        Raises:
            `ValueError`: When attempting to read a HMM from a closed
                file, or when the file could not be parsed.
            `~pyhmmer.errors.AllocationError`: When memory for the HMM could
                not be allocated successfully.
            `~pyhmmer.errors.AlphabetMismatch`: When the file contains HMMs
                in different alphabets, or in an alphabet that is different
                from the alphabet used to initialize the `HMMFile`.

        .. versionadded:: 0.4.11

        """
        cdef int   status
        cdef CM    py_cm
        cdef CM_t* cm     = NULL

        if self._fp == NULL:
            raise ValueError("I/O operation on closed file.")

        # don't run in *nogil* because the file may call a file-like handle
        status = libinfernal.cm_file.cm_file_Read(self._fp, True, &self._alphabet._abc, &cm)

        if status == libeasel.eslOK:
            return CM.from_ptr(cm, self._alphabet)
        elif status == libeasel.eslEOF:
            return None
        elif status == libeasel.eslEMEM:
            raise AllocationError("P7_HMM", sizeof(P7_HMM))
        elif status == libeasel.eslESYS:
            raise OSError(self._fp.errbuf.decode("utf-8", "replace"))
        elif status == libeasel.eslEFORMAT:
            raise ValueError("Invalid format in file: {}".format(self._fp.errbuf.decode("utf-8", "replace")))
        elif status == libeasel.eslEINCOMPAT:
            raise AlphabetMismatch(self._alphabet)
        else:
            _reraise_error()
            raise UnexpectedError(status, "p7_hmmfile_Read")

    cpdef void close(self) except *:
        """Close the CM file and free resources.

        This method has no effect if the file is already closed. It is called
        automatically if the `CMFile` was used in a context::

            >>> with CMFile("tests/data/cms/5.c.cm") as cm_file:
            ...     cm = cm_file.read()

        """
        if self._fp:
            libinfernal.cm_file.cm_file_Close(self._fp)
            self._fp = NULL


# cdef double   DEFAULT_F1      = 0.02
# cdef double   DEFAULT_F2      = 1e-3
# cdef double   DEFAULT_F3      = 1e-5
cdef uint32_t DEFAULT_SEED    = 181
# cdef double   DEFAULT_E       = 10.0
# cdef double   DEFAULT_INCE    = 0.01
# cdef size_t   HMMER_TARGET_LIMIT = 100000

cdef class Pipeline:

    CLEN_HINT = 100  # default model size
    L_HINT    = 100  # default sequence size

    cdef CM_PIPELINE* _pli
    cdef uint32_t     _seed
    cdef int64_t      _Z

    cdef readonly Alphabet   alphabet
    cdef readonly Randomness randomness

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._pli = NULL
        self.alphabet = None
        self.randomness = None

    def __init__(
        self,
        Alphabet alphabet,
        int64_t Z,
        *,
    #     bint bias_filter=True,
    #     bint null2=True,
        uint32_t seed=DEFAULT_SEED,
    #     object Z=None,
    #     object domZ=None,
    #     double F1=DEFAULT_F1,
    #     double F2=DEFAULT_F2,
    #     double F3=DEFAULT_F3,
    #     double E=DEFAULT_E,
    #     object T=None,
    #     double domE=DEFAULT_DOME,
    #     object domT=None,
    #     double incE=DEFAULT_INCE,
    #     object incT=None,
    #     double incdomE=DEFAULT_INCDOME,
    #     object incdomT=None,
    #     str bit_cutoffs=None,
    ):
        cdef int clen_hint = self.CLEN_HINT
        cdef int l_hint    = self.L_HINT

        with nogil:
            self._pli = libinfernal.cm_pipeline.cm_pipeline_Create(
                NULL,                               # ESL_GETOPTS *go
                alphabet._abc,                      # ESL_ALPHABET *abc
                clen_hint,                          # int clen_hint
                l_hint,                             # int L_hint
                Z,                                  # int Z
                cm_zsetby_e.CM_ZSETBY_OPTION,       # cm_zsetby_e Z_setby
                cm_pipemodes_e.CM_SEARCH_SEQS,      # cm_pipemodes_e mode
            )
        if self._pli == NULL:
            raise AllocationError("CM_PIPELINE", sizeof(CM_PIPELINE))

        # record alphabet
        self.alphabet = alphabet

        # create a Randomness object exposing the internal pipeline RNG
        self.randomness = Randomness.__new__(Randomness)
        self.randomness._owner = self
        self.randomness._rng = self._pli.r

        # configure the pipeline with the additional keyword arguments
        self.seed = seed
        self.Z = Z
        # self.E = E
        # self.T = T
        # self.incE = incE
        # self.incT = incT

        # self._pli.do_bot = False # FIXME: testing stuff

    def __dealloc__(self):
        # NOTE(@althonos): `cm_pipeline_Destroy` supposedly requires a `CM_t`
        #                  but does not use it so it *should* be fine to pass
        #                  a NULL pointer here.
        libinfernal.cm_pipeline.cm_pipeline_Destroy(self._pli, NULL)

    # --- Properties ---------------------------------------------------------

    @property
    def Z(self):
        """`int` or `None`: The number of effective targets searched.

        It is used to compute the independent e-value for each domain, and
        for an entire hit. If `None`, the parameter number will be set
        automatically after all the comparisons have been done. Otherwise,
        it can be set to an arbitrary number.

        """
        return None if self._Z < 0 else self._Z

    @Z.setter
    def Z(self, int64_t Z):
        assert self._pli != NULL
        if Z is None:
            self._pli.Z       = 0.0
            self._pli.Z_setby = cm_zsetby_e.CM_ZSETBY_OPTION
            self._Z           = -1
        else:
            self._pli.Z_setby = cm_zsetby_e.CM_ZSETBY_OPTION
            self._pli.Z = self._Z = Z

    @property
    def seed(self):
        """`int`: The seed given at pipeline initialization.

        Setting this attribute to a different value will cause the random
        number generator to be reseeded immediately.

        """
        return self._seed

    @seed.setter
    def seed(self, uint32_t seed):
        self._seed = seed
        self._pli.do_reseeding = self._pli.ddef.do_reseeding = seed != 0
        self.randomness.seed(seed)

    @property
    def F1(self):
        """`float`: The MSV filter threshold.
        """
        assert self._pli != NULL
        return self._pli.F1

    @F1.setter
    def F1(self, double F1):
        assert self._pli != NULL
        self._pli.F1 = F1

    @property
    def F2(self):
        """`float`: The Viterbi filter threshold.
        """
        assert self._pli != NULL
        return self._pli.F2

    @F2.setter
    def F2(self, double F2):
        assert self._pli != NULL
        self._pli.F2 = F2

    @property
    def F3(self):
        """`float`: The uncorrected Forward filter threshold.
        """
        assert self._pli != NULL
        return self._pli.F3

    @F3.setter
    def F3(self, double F3):
        assert self._pli != NULL
        self._pli.F3 = F3

    @property
    def F4(self):
        """`float`: The glocal Forward filter threshold.
        """
        assert self._pli != NULL
        return self._pli.F4

    @F4.setter
    def F4(self, double F4):
        assert self._pli != NULL
        self._pli.F4 = F4

    @property
    def F5(self):
        """`float`: The glocal envelope definition filter threshold.
        """
        assert self._pli != NULL
        return self._pli.F5

    @F5.setter
    def F5(self, double F5):
        assert self._pli != NULL
        self._pli.F5 = F5

    @property
    def F6(self):
        """`float`: The CYK filter threshold.
        """
        assert self._pli != NULL
        return self._pli.F6

    @F6.setter
    def F6(self, double F6):
        assert self._pli != NULL
        self._pli.F6 = F6

    @property
    def E(self):
        """`float`: The per-target E-value threshold for reporting a hit.
        """
        assert self._pli != NULL
        return self._pli.E

    @E.setter
    def E(self, double E):
        assert self._pli != NULL
        self._pli.E = E

    @property
    def T(self):
        """`float` or `None`: The per-target score threshold for reporting a hit.

        If set to a non-`None` value, this threshold takes precedence over
        the per-target E-value threshold (`Pipeline.E`).

        """
        assert self._pli != NULL
        return None if self._pli.by_E else self._pli.T

    @T.setter
    def T(self, object T):
        assert self._pli != NULL
        if T is None:
            self._pli.T = 0.0
            self._pli.by_E = True
        else:
            self._pli.T = T
            self._pli.by_E = False

    @property
    def incE(self):
        """`float`: The per-target E-value threshold for including a hit.

        .. versionadded:: 0.4.6

        """
        assert self._pli != NULL
        return self._pli.incE

    @incE.setter
    def incE(self, double incE):
        assert self._pli != NULL
        self._pli.incE = incE

    @property
    def incT(self):
        """`float` or `None`: The per-target score threshold for including a hit.

        If set to a non-`None` value, this threshold takes precedence over
        the per-target E-value inclusion threshold (`Pipeline.incE`).

        .. versionadded:: 0.4.8

        """
        assert self._pli != NULL
        return None if self._pli.inc_by_E else self._pli.incT

    @incT.setter
    def incT(self, object incT):
        assert self._pli != NULL
        if incT is None:
            self._pli.incT = 0.0
            self._pli.inc_by_E = True
        else:
            self._pli.incT = incT
            self._pli.inc_by_E = False

    # --- Utils --------------------------------------------------------------

    cdef int _configure_cm(
        self,
        WORKER_INFO* info,
    ) noexcept nogil:
        cdef int   status
        cdef float reqMb  = 0.0
        cdef bint check_fcyk_beta
        cdef bint check_final_beta
        cdef int  W_from_cmdline

        if (info.pli.cm_config_opts & libinfernal.cm.CM_CONFIG_SCANMX) != 0:
            reqMb += libinfernal.cm_mx.cm_scan_mx_SizeNeeded(info.cm, True, True)
        if (info.pli.cm_config_opts & libinfernal.cm.CM_CONFIG_TRSCANMX) != 0:
            reqMb += libinfernal.cm_mx.cm_tr_scan_mx_SizeNeeded(info.cm, True, True)
        if reqMb > info.smxsize:
            return libeasel.eslERANGE
            # ESL_FAIL(eslERANGE, info->pli->errbuf, "search will require %.2f Mb > %.2f Mb limit.\nIncrease limit with --smxsize, or don't use --max,--nohmm,--qdb,--fqdb.", reqMb, info->smxsize);

        # cm_pipeline_Create() sets configure/align options in pli->cm_config_opts, pli->cm_align_opts
        # NB: is this really necessary? we want to avoid modifying CM which is the query
        #     if we can help it
        info.cm.config_opts = info.pli.cm_config_opts
        info.cm.align_opts  = info.pli.cm_align_opts

        # check if we need to recalculate QDBs prior to building the scan matrix in cm_Configure()
        check_fcyk_beta  = (info.pli.fcyk_cm_search_opts & libinfernal.cm.CM_SEARCH_QDB) != 0
        check_final_beta = (info.pli.final_cm_search_opts & libinfernal.cm.CM_SEARCH_QDB) != 0
        if libinfernal.cm_qdband.CheckCMQDBInfo(info.cm.qdbinfo, info.pli.fcyk_beta, check_fcyk_beta, info.pli.final_beta, check_final_beta) != libeasel.eslOK:
            info.cm.config_opts  |= libinfernal.cm.CM_CONFIG_QDB
            info.cm.qdbinfo.beta1 = info.pli.fcyk_beta
            info.cm.qdbinfo.beta2 = info.pli.final_beta

        W_from_cmdline = -1 if not info.pli.do_wcx else <int> (info.cm.clen * info.pli.wcx)
        return libinfernal.cm_modelconfig.cm_Configure(info.cm, info.pli.errbuf, W_from_cmdline)

    cdef int _setup_hmm_filter(
        self,
        WORKER_INFO* info,
    ) except * nogil:
        cdef bint do_trunc_ends = True #(esl_opt_GetBoolean(go, "--notrunc") || esl_opt_GetBoolean(go, "--inttrunc")) ? FALSE : TRUE;

        # set up the HMM filter-related structures
        info.gm = libhmmer.p7_profile.p7_profile_Create(info.cm.fp7.M, info.cm.abc)
        info.om = p7_oprofile_Create(info.cm.fp7.M, info.cm.abc)
        info.bg = libhmmer.p7_bg.p7_bg_Create(info.cm.abc)
        libhmmer.modelconfig.p7_ProfileConfig(info.cm.fp7, info.bg, info.gm, 100, libhmmer.p7_LOCAL) # 100 is a dummy length for now; and MSVFilter requires local mode
        p7_oprofile_Convert(info.gm, info.om)                          # <om> is now p7_LOCAL, multihit
        # clone gm into Tgm before putting it into glocal mode
        if do_trunc_ends:
            info.Tgm = libhmmer.p7_profile.p7_profile_Clone(info.gm)

        # after om has been created, convert gm to glocal, to define envelopes in cm_pipeline()
        libhmmer.modelconfig.p7_ProfileConfig(info.cm.fp7, info.bg, info.gm, 100, libhmmer.p7_GLOCAL)

        if do_trunc_ends:
            # create Rgm, Lgm, and Tgm specially-configured profiles for defining envelopes around
            # hits that may be truncated 5' (Rgm), 3' (Lgm) or both (Tgm).
            info.Rgm = libhmmer.p7_profile.p7_profile_Clone(info.gm)
            info.Lgm = libhmmer.p7_profile.p7_profile_Clone(info.gm)
            # info->Tgm was created when gm was still in local mode above
            # we cloned Tgm from the while profile was still locally configured, above
            libinfernal.cm_p7_modelconfig.p7_ProfileConfig5PrimeTrunc(info.Rgm, 100)
            libinfernal.cm_p7_modelconfig.p7_ProfileConfig3PrimeTrunc(info.cm.fp7, info.Lgm, 100)
            libinfernal.cm_p7_modelconfig.p7_ProfileConfig5PrimeAnd3PrimeTrunc(info.Tgm, 100)
        else:
            info.Rgm = NULL
            info.Lgm = NULL
            info.Tgm = NULL

        # copy E-value parameters
        libeasel.vec.esl_vec_FCopy(info.cm.fp7_evparam, libinfernal.CM_p7_NEVPARAM, info.p7_evparam);

        # compute msvdata
        info.msvdata = libhmmer.p7_scoredata.p7_hmm_ScoreDataCreate(info.om, NULL)

        return libeasel.eslOK

    cdef void free_info(
        self,
        WORKER_INFO *info,
    ) noexcept nogil:
        # TODO: free or use Python garbage collection with dedicated objects
        # if(info.pli != NULL):
        #     cm_pipeline_Destroy(info.pli, info.cm);
        #     info.pli = NULL
        # if(info.th != NULL):
        #     cm_tophits_Destroy(info->th)
        #     info->th = NULL
        # if(info.cm != NULL):
        #     FreeCM(info->cm)
        #     info->cm = NULL

        # if(info.om != NULL):
        #     p7_oprofile_Destroy(info->om)
        #     info->om = NULL
        # if(info.gm != NULL):
        #     p7_profile_Destroy(info->gm)
        #     info->gm = NULL
        # if(info.Rgm != NULL):
        #     p7_profile_Destroy(info->Rgm)
        #     info->Rgm = NULL
        # if(info.Lgm != NULL):
        #     p7_profile_Destroy(info->Lgm);
        #     info->Lgm = NULL
        # if(info.Tgm != NULL):
        #     p7_profile_Destroy(info->Tgm);
        #     info->Tgm = NULL
        # if(info.bg != NULL):
        #     p7_bg_Destroy(info->bg)
        #     info->bg = NULL
        # if(info.p7_evparam != NULL):
        #     free(info->p7_evparam)
        #     info->p7_evparam = NULL
        # if(info.msvdata    != NULL):
        #     p7_hmm_ScoreDataDestroy(info->msvdata);
        #     info->msvdata    = NULL
        return


    # --- Methods ------------------------------------------------------------

    @staticmethod
    cdef int _search_loop(
        WORKER_INFO* info,
        ESL_SQ** sq,
        size_t n_targets,
        int nbps,
    ) except 1 nogil:
        # adapted from `serial_loop` in `cmsearch.c`, inner loop code

        cdef int      status
        cdef size_t   t
        cdef uint64_t prv_pli_ntophits

        # prepare pipeline for new model
        status = libinfernal.cm_pipeline.cm_pli_NewModel(
            info.pli,
            cm_newmodelmodes_e.CM_NEWMODEL_CM,
            info.cm,
            info.cm.clen,
            info.cm.W,
            nbps,
            info.om,
            info.bg,
            info.p7_evparam,
            info.om.max_length,
            0, #cm_idx - 1, # FIXME?  # int64_t cur_cm_idx
            -1,                       # int     cur_clan_idx
            NULL,
        )
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "cm_pli_NewModel")

        # run the inner loop on all sequences
        for t in range(n_targets):
            # configure the pipeline for a new sequence
            status = libinfernal.cm_pipeline.cm_pli_NewSeq(info.pli, sq[t], t)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "cm_pli_NewSeq")

            # run top strand
            if info.pli.do_top:
                prv_pli_ntophits = info.th.N
                status = libinfernal.cm_pipeline.cm_Pipeline(info.pli, info.cm.offset, info.om, info.bg, info.p7_evparam, info.msvdata, sq[t], info.th, False, NULL, &info.gm, &info.Rgm, &info.Lgm, &info.Tgm, &info.cm)
                if status != libeasel.eslOK:
                    raise EaselError(status, info.pli.errbuf.decode('utf-8', 'ignore'))
                libinfernal.cm_pipeline.cm_pipeline_Reuse(info.pli)  # prepare for next search
                if sq[t].C > 0:
                    libinfernal.cm_pipeline.cm_pli_AdjustNresForOverlaps(info.pli, sq[t].C, False)
                libinfernal.cm_tophits.cm_tophits_UpdateHitPositions(info.th, prv_pli_ntophits, sq[t].start, False)

            # reverse complement
            if info.pli.do_bot and sq[t].abc.complement != NULL:
                prv_pli_ntophits = info.th.N
                libeasel.sq.esl_sq_ReverseComplement(sq[t])  # FIXME: make a copy?
                status = libinfernal.cm_pipeline.cm_Pipeline(info.pli, info.cm.offset, info.om, info.bg, info.p7_evparam, info.msvdata, sq[t], info.th, True, NULL, &info.gm, &info.Rgm, &info.Lgm, &info.Tgm, &info.cm)
                if status != libeasel.eslOK:
                    raise EaselError(status, info.pli.errbuf.decode('utf-8', 'ignore'))
                libinfernal.cm_pipeline.cm_pipeline_Reuse(info.pli)  # prepare for next search
                if sq[t].C > 0:
                    libinfernal.cm_pipeline.cm_pli_AdjustNresForOverlaps(info.pli, sq[t].C, True)
                libinfernal.cm_tophits.cm_tophits_UpdateHitPositions(info.th, prv_pli_ntophits, sq[t].start, True)
                libeasel.sq.esl_sq_ReverseComplement(sq[t])

        # Return 0 to indicate success
        return 0

    cpdef TopHits search_cm(
        self,
        CM query,
        SearchTargets sequences,
    ):
        # adapted from `serial_master` in `cmsearch.c`, outer loop code

        cdef float[CM_p7_NEVPARAM] p7_evparam
        cdef WORKER_INFO           tinfo
        cdef int                   status
        cdef double                eZ
        cdef TopHits               top_hits = TopHits(query)

        # check that all alphabets are consistent
        if not self.alphabet._eq(query.alphabet):
            raise AlphabetMismatch(self.alphabet, query.alphabet)
        if not self.alphabet._eq(sequences.alphabet):
            raise AlphabetMismatch(self.alphabet, sequences.alphabet)

        tinfo.p7_evparam = p7_evparam
        tinfo.smxsize = 128.0
        tinfo.cm = query._cm
        tinfo.pli = self._pli  # Maybe copy?
        tinfo.th = top_hits._th

        if (tinfo.cm.flags & libinfernal.cm.CMH_FP7) == 0:
            raise ValueError(f"no filter HMM was found for CM {query.name!r}")

        # check if we have E-value stats for the CM, we require them
        # *unless* we are going to run the pipeline in HMM-only mode.
        # We run the pipeline in HMM-only mode if --nohmmonly is
        # not used and -g is not used and:
        # (a) --hmmonly used OR
        # (b) model has 0 basepairs
        nbps = libinfernal.cm.CMCountNodetype(tinfo.cm, libinfernal.MATP_nd)
        # TODO: below
        # if((   esl_opt_GetBoolean(go, "--nohmmonly"))  ||
        # (   esl_opt_GetBoolean(go, "-g"))           ||
        # ((! esl_opt_GetBoolean(go, "--hmmonly"))    && (nbps > 0))) {
        # /* we're NOT running HMM-only pipeline variant, we need CM E-value stats */
        # if(! (tinfo->cm->flags & CMH_EXPTAIL_STATS)) cm_Fail("no E-value parameters were read for CM: %s.\nYou may need to run cmcalibrate.", tinfo->cm->name);
        # }

        # configure the CM (this builds QDBs if nec) and setup HMM filters
        # (we need to do this before clone_info()). We need a pipeline to
        # do this only b/c we need pli->cm_config_opts.
        #
        status = self._configure_cm(&tinfo)
        if status != libeasel.eslOK:
            raise EaselError(status, tinfo.pli.errbuf.decode('utf-8', 'ignore'))
        status = self._setup_hmm_filter(&tinfo)
        if status != libeasel.eslOK:
            raise EaselError(status, tinfo.pli.errbuf.decode('utf-8', 'ignore'))

        with nogil:
            # make sure the pipeline is set to search mode
            self._pli.mode = cm_pipemodes_e.CM_SEARCH_SEQS
            # run the cmsearch loop on all database sequences while
            # recycling memory between targets
            if SearchTargets is DigitalSequenceBlock:
                Pipeline._search_loop(&tinfo, sequences._refs, sequences._length, nbps)
            elif SearchTargets is SequenceFile:
                raise NotImplementedError("Pipeline.search_cm")
            else:
                raise NotImplementedError("Pipeline.search_cm")

        # we need to re-compute e-values before merging (when list will be sorted)
        if tinfo.pli.do_hmmonly_cur:
            eZ = tinfo.pli.Z / <float> tinfo.om.max_length
        else:
            eZ = tinfo.cm.expA[tinfo.pli.final_cm_exp_mode].cur_eff_dbsize
        libinfernal.cm_tophits.cm_tophits_ComputeEvalues(tinfo.th, eZ, 0)

        ## merge the search results
        # cm_tophits_Merge(t.th,   info.th);
        # cm_pipeline_Merge(info[0].pli, info.pli);
        # free_info(&info)

        # Sort by sequence index/position and remove duplicates
        libinfernal.cm_tophits.cm_tophits_SortForOverlapRemoval(tinfo.th)
        status = libinfernal.cm_tophits.cm_tophits_RemoveOrMarkOverlaps(tinfo.th, False, tinfo.pli.errbuf)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "cm_tophits_RemoveOrMarkOverlaps")

        # Resort: by score (usually) or by position (if in special 'terminate after F3' mode) */
        if tinfo.pli.do_trm_F3:
            libinfernal.cm_tophits.cm_tophits_SortByPosition(tinfo.th)
        else:
            libinfernal.cm_tophits.cm_tophits_SortByEvalue(tinfo.th)

        # Enforce threshold
        libinfernal.cm_tophits.cm_tophits_Threshold(tinfo.th, tinfo.pli)

        # tally up total number of hits and target coverage
        # NOTE: likely uneeded as we don't report accounting?
        # for i in range(tinfo.th.N):
        #     if (tinfo.th.hit[i].flags & (libinfernal.cm_tophits.CM_HIT_IS_REPORTED | libinfernal.cm_tophits.CM_HIT_IS_REPORTED)) != 0:
        #         tinfo.pli.acct[tinfo.th.hit[i].pass_idx].n_output += 1
        #         tinfo.pli.acct[tinfo.th.hit[i].pass_idx].pos_output += abs(tinfo.th.hit[i].stop - tinfo.th.hit[i].start) + 1

        # Record pipeline configuration before returning
        memcpy(&top_hits._pli, self._pli, sizeof(CM_PIPELINE)) # FIXME?
        return top_hits


cdef class Hit:
    # a reference to the TopHits that owns the wrapped CM_HIT, kept so that
    # the internal data is never deallocated before the Python class.
    cdef readonly TopHits hits
    cdef CM_HIT* _hit

    def __cinit__(self, TopHits hits, size_t index):
        assert hits._th != NULL
        assert index < hits._th.N
        self.hits = hits
        self._hit = hits._th.hit[index]

    @property
    def name(self):
        """`bytes`: The name of the database hit.
        """
        assert self._hit != NULL
        assert self._hit.name != NULL
        return <bytes> self._hit.name

    @name.setter
    def name(self, bytes name not None):
        assert self._hit != NULL
        free(self._hit.name)
        self._hit.name = strdup(<const char*> name)
        if self._hit.name == NULL:
            raise AllocationError("char", sizeof(char), strlen(name))


cdef class TopHits:
    cdef CM_TOPHITS* _th
    cdef CM_PIPELINE _pli
    cdef object      _query

    def __cinit__(self):
        self._th = NULL
        self._query = None
        memset(&self._pli, 0, sizeof(CM_PIPELINE))

    def __init__(self, object query not None):
        self._query = query
        with nogil:
            # free allocated memory (in case __init__ is called more than once)
            libinfernal.cm_tophits.cm_tophits_Destroy(self._th)
            # allocate top hits
            self._th = libinfernal.cm_tophits.cm_tophits_Create()
            if self._th == NULL:
                raise AllocationError("CM_TOPHITS", sizeof(CM_TOPHITS))
            # clear pipeline configuration
            memset(&self._pli, 0, sizeof(CM_PIPELINE))

    def __dealloc__(self):
        libinfernal.cm_tophits.cm_tophits_Destroy(self._th)

    def __bool__(self):
        assert self._th != NULL
        return self._th.N > 0

    def __len__(self):
        assert self._th != NULL
        return self._th.N

    def __getitem__(self, index):
        assert self._th != NULL
        if not (
               self._th.is_sorted_by_evalue
            or self._th.is_sorted_for_overlap_removal
            or self._th.is_sorted_for_overlap_markup
            or self._th.is_sorted_by_position
        ):
            for i in range(self._th.N):
                self._th.hit[i] = &self._th.unsrt[i]
        if index < 0:
            index += self._th.N
        if index >= self._th.N or index < 0:
            raise IndexError("list index out of range")
        return Hit(self, index)

    # --- Properties ---------------------------------------------------------

    @property
    def query(self):
        """`object`: The query object these hits were obtained for.

        The actual type of `TopHits.query` depends on the query that was given
        to the `Pipeline` that created the object.

        """
        return self._query

    # --- Methods ------------------------------------------------------------

    cpdef void write(self, object fh, str format="3", bint header=True) except *:
        """Write the hits in tabular format to a file-like object.

        Arguments:
            fh (`io.IOBase`): A Python file handle, opened in binary mode.
            format (`str`): The tabular format in which to write the hits.
            header (`bool`): Whether to write a table header. Ignored
                when writing in the ``pfam`` format.

        """
        cdef FILE* file
        cdef str   fname
        cdef int   status
        cdef bytes qname  = b"-"
        cdef bytes qacc   = b"-"

        if self._query is not None:
            if self._query.name is not None:
                qname = self._query.name
            if self._query.accession is not None:
                qacc = self._query.accession

        file = fopen_obj(fh, "w")
        try:
            if format == "3":
                fname = "cm_tophits_TabularTargets3"
                status = libinfernal.cm_tophits.cm_tophits_TabularTargets3(
                    file,
                    qname,
                    qacc,
                    self._th,
                    &self._pli,
                    header
                )
            else:
                raise InvalidParameter("format", format, choices=["3"])
            if status != libeasel.eslOK:
                _reraise_error()
                raise UnexpectedError(status, fname)
        finally:
            fclose(file)

# --- Module init code -------------------------------------------------------

init_ilogsum()
FLogsumInit()
p7_FLogsumInit()