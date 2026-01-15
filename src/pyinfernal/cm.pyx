from cpython.exc cimport PyErr_WarnEx
from cpython.unicode cimport (
    PyUnicode_FromString,
    PyUnicode_DecodeASCII,
    PyUnicode_AsUTF8AndSize,
)

from libc cimport errno
from libc.stdio cimport FILE, fopen, fclose
from libc.stdint cimport uint32_t, uint64_t, int64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memset, memcpy, memmove, strdup, strndup, strncpy, strlen

cimport libeasel
cimport libeasel.alphabet
cimport libeasel.vec
cimport libeasel.fileparser
cimport libhmmer.p7_bg
cimport libhmmer.p7_profile
cimport libhmmer.p7_hmm
cimport libhmmer.p7_gmx
cimport libhmmer.p7_domaindef
cimport libhmmer.p7_scoredata
cimport libhmmer.modelconfig
cimport libinfernal.cm
cimport libinfernal.cm_alidisplay
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
from libhmmer.p7_scoredata cimport P7_SCOREDATA
from libhmmer.p7_gmx cimport P7_GMX
from libhmmer.logsum cimport p7_FLogsumInit
from libinfernal cimport CM_p7_NEVPARAM
from libinfernal.cm_file cimport CM_FILE, cm_file_formats_e
from libinfernal.cm_pipeline cimport CM_PIPELINE, CM_PLI_ACCT, cm_zsetby_e, cm_pipemodes_e, cm_newmodelmodes_e
from libinfernal.cm_tophits cimport CM_TOPHITS, CM_HIT
from libinfernal.cm cimport CM_t
from libinfernal.cmsearch cimport WORKER_INFO
from libinfernal.logsum cimport FLogsumInit, init_ilogsum
from libinfernal.cm_alidisplay cimport CM_ALIDISPLAY

if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_omx cimport P7_OM_BLOCK, p7_omx_Reuse
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE, p7_oprofile_Create, p7_oprofile_Convert
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_omx cimport P7_OM_BLOCK, p7_omx_Reuse
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE, p7_oprofile_Create, p7_oprofile_Convert
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon.p7_omx cimport P7_OM_BLOCK, p7_omx_Reuse
    from libhmmer.impl_neon.p7_oprofile cimport P7_OPROFILE, p7_oprofile_Create, p7_oprofile_Convert

if TARGET_SYSTEM == "Linux":
    from .fileobj.linux cimport fileobj_linux_open as fopen_obj
elif TARGET_SYSTEM == "Darwin" or TARGET_SYSTEM.endswith("BSD"):
    from .fileobj.bsd cimport fileobj_bsd_open as fopen_obj

cimport pyhmmer.easel
cimport pyhmmer.plan7
from pyhmmer.easel cimport (
    Alphabet,
    DigitalSequenceBlock,
    Randomness,
    SequenceFile,
)
from pyhmmer.plan7 cimport (
    Background,
    HMM,
    Profile,
    OptimizedProfile,
)

include "exceptions.pxi"
include "_strings.pxi"
# include "_getid.pxi"

# --- Python imports ---------------------------------------------------------

import datetime
import enum
import io
import os
import operator
import sys
import warnings

from pyhmmer.utils import SizedIterator
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

    Attributes:
        filter_hmm (`pyhmmer.plan7.HMM` or `None`): The HMM used for
            the initial filtering stages inside the Infernal pipeline.

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
            obj.alphabet = Alphabet.from_ptr(cm.abc)
        else:
            obj.alphabet = alphabet

        if cm.flags & libinfernal.cm.CMH_FP7:
            obj.filter_hmm = HMM.__new__(HMM)
            obj.filter_hmm._hmm = cm.fp7
            obj.filter_hmm.alphabet = obj.alphabet

        if cm.flags & libinfernal.cm.CMH_MLP7:
            obj.ml_hmm = HMM.__new__(HMM)
            obj.ml_hmm._hmm = cm.mlp7
            obj.ml_hmm.alphabet = obj.alphabet

        return obj

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self.alphabet = None
        self.filter_hmm = None
        self.ml_hmm = None
        self._cm = NULL

    def __init__(
        self,
        Alphabet alphabet not None,
        int M,
        int N,
        int clen,
        str name not None,
    ):
        """__init__(self, alphabet, M, name)\n--\n

        Create a new, empty CM from scratch.

        Arguments:
            alphabet (`~pyhmmer.easel.Alphabet`): The alphabet of the model.
            M (`int`): The number of nodes of the model.
            N (`int`): The number of states of the model.
            clen (`int`): The model consensus length.
            name (`str`): The name of the model.

        """
        # store the alphabet so it's not deallocated
        self.alphabet = alphabet
        # create a new HMM suitable for at least M nodes
        with nogil:
            self._cm = libinfernal.cm.CreateCM(M, N, clen, alphabet._abc)
        if not self._cm:
            raise AllocationError("CM_t", sizeof(CM_t))
        # record mandatory name
        self.name = name

    def __dealloc__(self):
        if self._cm is not NULL:
            self._cm.fp7 = NULL # owned by `self.filter_hmm`
            self._cm.mlp7 = NULL # owned by `self.ml_hmm`
            libinfernal.cm.FreeCM(self._cm)

    def __sizeof__(self):
        assert self._cm != NULL
        return sizeof(self) + libinfernal.cm.cm_Sizeof(self._cm)

    # --- Properties ---------------------------------------------------------

    @property
    def N(self):
        """`int`: The number of nodes in the model.
        """
        assert self._cm != NULL
        return self._cm.nodes

    @property
    def M(self):
        """`int`: The number of states in the model.
        """
        assert self._cm != NULL
        return self._cm.M

    @property
    def name(self):
        """`str`: The name of the CM.
        """
        assert self._cm != NULL
        return PyUnicode_FromString(self._cm.name)

    @name.setter
    def name(self, object name not None):
        assert self._cm != NULL
        status = _set_str(self._cm, name, <str_setter_t> libinfernal.cm.cm_SetName)
        if status == libeasel.eslEMEM:
            raise AllocationError("char", sizeof(char), len(name))
        elif status != libeasel.eslOK:
            raise UnexpectedError(status, "cm_SetName")

    @property
    def accession(self):
        """`str` or `None`: The accession of the CM, if any.
        """
        assert self._cm != NULL
        if not self._cm.flags & libinfernal.cm.CMH_ACC:
            return None
        assert self._cm.acc != NULL
        return PyUnicode_FromString(self._cm.acc)

    @property
    def description(self):
        """`str` or `None`: The description of the CM, if any.
        """
        assert self._cm != NULL
        if not self._cm.flags & libinfernal.cm.CMH_DESC:
            return None
        assert self._cm.desc != NULL
        return PyUnicode_FromString(self._cm.desc)

    @property
    def command_line(self):
        """`str` or `None`: The command line that built the model.
        """
        assert self._cm != NULL
        if self._cm.comlog == NULL:
            return None
        return self._cm.comlog.decode("ascii")

    @command_line.setter
    def command_line(self, str cli):
        assert self._cm != NULL

        cdef const char*  cli_
        cdef ssize_t      n    = -1

        if cli is None:
            free(self._cm.comlog)
            self._cm.comlog = NULL
        else:
            cli_ = PyUnicode_AsUTF8AndSize(cli, &n)
            self._cm.comlog = <char*> realloc(<void*> self._cm.comlog, sizeof(char) * (n + 1))
            if self._cm.comlog == NULL:
                raise AllocationError("char", sizeof(char), n+1)
            with nogil:
                strncpy(self._cm.comlog, cli_, n+1)

    @property
    def nseq(self):
        """`int` or `None`: The number of training sequences used, if any.
        """
        assert self._cm != NULL
        return None if self._cm.nseq == -1 else self._cm.nseq

    @property
    def nseq_effective(self):
        """`float` or `None`: The number of effective sequences used, if any.
        """
        assert self._cm != NULL
        return None if self._cm.eff_nseq == -1.0 else self._cm.eff_nseq

    @property
    def creation_time(self):
        """`datetime.datetime` or `None`: The creation time of the CM, if any.

        Example:
            Get the creation time for any HMM::

                >>> lsu_rrna.creation_time
                datetime.datetime(2014, 6, 16, 10, 10)

            Set the creation time manually to a different date and time::

                >>> ctime = datetime.datetime(2026, 1, 15, 13, 59, 10)
                >>> lsu_rrna.creation_time = ctime
                >>> lsu_rrna.creation_time
                datetime.datetime(2026, 1, 15, 13, 59, 10)

        Danger:
            Internally, Infernal always uses ``asctime`` to generate a
            timestamp for the HMMs, so this property assumes that every
            creation time field can be parsed into a `datetime.datetime`
            object using the  ``"%a %b %d %H:%M:%S %Y"`` format.

        """
        assert self._cm != NULL

        cdef size_t l
        cdef str    ctime

        if self._cm.ctime == NULL:
            return None

        l = strlen(self._cm.ctime)
        ctime = PyUnicode_DecodeASCII(self._cm.ctime, l, NULL)
        return datetime.datetime.strptime(ctime,'%a %b %d %H:%M:%S %Y')

    @creation_time.setter
    def creation_time(self, object ctime):
        assert self._cm != NULL

        cdef str         ty
        cdef str         formatted
        cdef const char* s
        cdef ssize_t     n

        if ctime is None:
            free(self._cm.ctime)
            self._cm.ctime = NULL
            return
        elif not isinstance(ctime, datetime.datetime):
            ty = type(ctime).__name__
            raise TypeError(f"Expected datetime.datetime or None, found {ty}")

        formatted = ctime.strftime('%a %b %e %H:%M:%S %Y')
        s = PyUnicode_AsUTF8AndSize(formatted, &n)

        self._cm.ctime = <char*> realloc(<void*> self._cm.ctime, sizeof(char) * (n + 1))
        if self._cm.ctime == NULL:
            raise AllocationError("char", sizeof(char), n+1)
        with nogil:
            strncpy(self._cm.ctime, s, n + 1)


cdef class CMFile:
    """A wrapper around a file storing serialized CMs.

    Example:
        Load the first CM from a CM file located on the
        local filesystem::

            >>> with CMFile("tests/data/cms/iss33.cm") as cm_file:
            ...     cm = cm_file.read()
            >>> cm.name
            'LSU_rRNA_bacteria'
            >>> cm.accession
            'RF02541'

        Load all the CMs from a CM file into a `list`::

            >>> with CMFile("tests/data/cms/5.c.cm") as cm_file:
            ...     cms = list(cm_file)
            >>> len(cms)
            5
            >>> [cm.accession for cm in cms]
            ['RF00005', 'RF00006', 'RF01185', 'RF01855', 'RF01852']

    """

    cdef CM_FILE* _fp
    cdef str      _name
    cdef Alphabet _alphabet
    cdef ESL_ALPHABET* _abc

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
        self._abc = NULL
        self._fp = NULL
        self._name = None

    def __init__(self, object file, bint db = True, *, Alphabet alphabet = None):
        """__init__(self, file, db=True, *, alphabet=None)\n--\n

        Create a CM reader from the given path or file.

        Arguments:
            file (`str`, `bytes`, `os.PathLike` or file-like object): Either
                the path to a file containing the CMs to read, or a file-like
                object in **binary mode**.
            db (`bool`): Set to `False` to force the parser to ignore the
                pressed HMM database if it finds one. Defaults to `True`.
            alphabet (`~pyhmmer.easel.Alphabet`, optional): The alphabet
                of the CMs in the file. Supports auto-detection, but passing
                a non-`None` argument will facilitate MyPy type inference.

        Raises:
            `TypeError`: When ``file`` is not of the correct type, or when
                ``file`` is a file-like object open in text mode rather
                than binary mode.
            `RuntimeError`: When the internal system function
                (``fopencookie`` on Linux, ``funopen`` on BSD) fails to open
                the file.

        """
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

        if alphabet is None:
            self._alphabet = None
            self._abc = NULL
        else:
            self._alphabet = alphabet
            self._abc = alphabet._abc

    def __dealloc__(self):
        if self._fp:
            PyErr_WarnEx(ResourceWarning, "unclosed CM file", 1)
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

        """
        cdef int   status
        cdef CM    py_cm
        cdef CM_t* cm     = NULL

        if self._fp == NULL:
            raise ValueError("I/O operation on closed file.")

        # don't run in *nogil* because the file may call a file-like handle
        status = libinfernal.cm_file.cm_file_Read(self._fp, True, &self._abc, &cm)

        if self._alphabet is None and self._abc != NULL:
            self._alphabet = Alphabet.from_ptr(self._abc)

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
            >>> cm_file.closed
            True

        """
        if self._fp:
            libinfernal.cm_file.cm_file_Close(self._fp)
            self._fp = NULL


cdef uint32_t DEFAULT_SEED    = 181
cdef double   DEFAULT_E       = 10.0
cdef double   DEFAULT_INCE    = 0.01

cdef class Pipeline:
    """An Infernal accelerated sequence/covariance model comparison pipeline.

    The Infernal pipeline handles the comparison of digital sequences to
    CMs. Since Infernal 1.1, the pipeline is accelerated using the HMMER
    platform acceleration to compute initial SSV/MSV/Viterbi filters before
    actually attempting full alignments.

    Attributes:
        alphabet (`~pyhmmer.easel.Alphabet`): The alphabet for which the
            pipeline is configured.
        randomness (`~pyhmmer.easel.Randomness`): The random number generator
            being used by the pipeline.

    """

    CLEN_HINT = 100  # default model size
    M_HINT    = 100  # default HMM size
    L_HINT    = 100  # default sequence size

    cdef CM_PIPELINE* _pli
    cdef uint32_t     _seed
    cdef int64_t      _Z

    cdef readonly Alphabet         alphabet
    cdef readonly Randomness       randomness
    cdef readonly Background       background

    cdef          OptimizedProfile opt          # temporary optimized profile
    cdef          Profile          profile      # temporary profile
    cdef          Profile          profile_l    # temporary profile, 3' truncated
    cdef          Profile          profile_r    # temporary profile, 5' truncated
    cdef          Profile          profile_t    # temporary profile, 5' + 3' truncated

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._pli = NULL
        self.alphabet = None
        self.randomness = None

    def __init__(
        self,
        Alphabet alphabet,
        int64_t Z,
        Background background = None,
        *,
        # bint bias_filter=True,
        # bint null2=True,
        uint32_t seed=DEFAULT_SEED,
        # object Z=None,
        # double F1=DEFAULT_F1,
        # double F2=DEFAULT_F2,
        # double F3=DEFAULT_F3,
        double E=DEFAULT_E,
        object T=None,
        double incE=DEFAULT_INCE,
        object incT=None,
    #     str bit_cutoffs=None,
    ):
        cdef int clen_hint = self.CLEN_HINT
        cdef int l_hint    = self.L_HINT
        cdef int m_hint    = self.M_HINT

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

        # use the backgroud model or create a default one
        if background is None:
            self.background = Background(alphabet)
        elif background.alphabet != self.alphabet:
            raise AlphabetMismatch(self.alphabet, background.alphabet)
        else:
            self.background = background.copy() # FIXME: do we really need a copy?

        # create a Randomness object exposing the internal pipeline RNG
        self.randomness = Randomness.__new__(Randomness)
        self.randomness._owner = self
        self.randomness._rng = self._pli.r

        # create empty profiles and optimized profile to reuse globally
        # between queries rather than reallocating on every new query
        self.profile = Profile(m_hint, self.alphabet)
        self.opt = OptimizedProfile(m_hint, self.alphabet)
        self.profile_r = Profile(m_hint, self.alphabet)
        self.profile_l = Profile(m_hint, self.alphabet)
        self.profile_t = Profile(m_hint, self.alphabet)

        # configure the pipeline with the additional keyword arguments
        self.seed = seed
        self.Z = Z
        self.E = E
        self.T = T
        self.incE = incE
        self.incT = incT

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

        It is used to compute the independent e-value for each hit.
        If `None`, the parameter number will be set automatically after all
        the comparisons have been done. Otherwise, it can be set to an
        arbitrary number.

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

    cpdef void clear(self):
        """Reset the pipeline to its default state.
        """
        assert self._pli != NULL

        cdef int      status
        cdef int      i
        cdef uint32_t seed

        # reinitialize the random number generator, even if
        # `self._pli.do_reseeding` is False, because a true
        # deallocation/reallocation of a P7_PIPELINE would reinitialize
        # it unconditionally.
        self.randomness.seed(self._seed)

        # reinitialize the domaindef
        status = libhmmer.p7_domaindef.p7_domaindef_Reuse(self._pli.ddef)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "p7_domaindef_Reuse")

        with nogil:
            # reinitialize the dynamic programming matrices
            # (no status check, infallible)
            p7_omx_Reuse(self._pli.oxf)
            p7_omx_Reuse(self._pli.oxb)
            p7_omx_Reuse(self._pli.fwd)
            p7_omx_Reuse(self._pli.bck)
            libhmmer.p7_gmx.p7_gmx_Reuse(self._pli.gxf)
            libhmmer.p7_gmx.p7_gmx_Reuse(self._pli.gxb)
            libhmmer.p7_gmx.p7_gmx_Reuse(self._pli.gfwd)
            libhmmer.p7_gmx.p7_gmx_Reuse(self._pli.gbck)
            # Reset accounting values
            self._pli.nseqs           = 0
            self._pli.nmodels         = 0
            self._pli.nnodes          = 0
            self._pli.nmodels_hmmonly = 0
            self._pli.nnodes_hmmonly  = 0
            self._pli.cmfp            = NULL
            self._pli.errbuf[0]       = b'\0'
            # Mass reset of pipeline accounting statistics
            # (we could use `cm_pli_ZeroAccounting`, but because every field
            # is a `uint64_t`, we can just set the whole array to zero bytes)
            memset(self._pli.acct, 0, libinfernal.cm_pipeline.NPLI_PASSES * sizeof(CM_PLI_ACCT))

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
        # NOTE: is this really necessary? we want to avoid modifying CM which is the query
        #       if we can help it
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

    cdef int _grow_profiles(
        self,
        CM query,
    ) except 1:
        cdef int M = query.M
        if self.profile._gm.allocM < M:
            self.profile = Profile(M, self.alphabet)
            self.opt = OptimizedProfile(M, self.alphabet)
            self.profile_r = Profile(M, self.alphabet)
            self.profile_l = Profile(M, self.alphabet)
            self.profile_t = Profile(M, self.alphabet)
        else:
            self.profile.clear()
            self.profile_r.clear()
            self.profile_l.clear()
            self.profile_t.clear()

    cdef int _setup_hmm_filter(
        self,
        WORKER_INFO* info,
        CM query,
    ) except 1:
        cdef int status
        cdef bint do_trunc_ends = True #(esl_opt_GetBoolean(go, "--notrunc") || esl_opt_GetBoolean(go, "--inttrunc")) ? FALSE : TRUE;

        # set up the HMM filter-related structures
        self._grow_profiles(query)
        info.gm = self.profile._gm
        info.om = self.opt._om
        info.bg = self.background._bg
        self.profile.configure(query.filter_hmm, self.background, 100, multihit=True, local=True)
        self.opt.convert(self.profile) # <om> is now p7_LOCAL, multihit

        # clone gm into Tgm before putting it into glocal mode
        if do_trunc_ends:
            status = libhmmer.p7_profile.p7_profile_Copy(self.profile._gm, self.profile_t._gm)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_profile_Copy")

        # after om has been created, convert gm to glocal, to define envelopes in cm_pipeline()
        self.profile.configure(query.filter_hmm, self.background, 100, multihit=True, local=False)

        if do_trunc_ends:
            # create Rgm, Lgm, and Tgm specially-configured profiles for defining envelopes around
            # hits that may be truncated 5' (Rgm), 3' (Lgm) or both (Tgm).
            # FIXME: reuse instead of copying (which causes reallocation)
            status = libhmmer.p7_profile.p7_profile_Copy(self.profile._gm, self.profile_r._gm)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_profile_Copy")
            status = libhmmer.p7_profile.p7_profile_Copy(self.profile._gm, self.profile_l._gm)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_profile_Copy")
            # setup pointers
            info.Tgm = self.profile_t._gm
            info.Rgm = self.profile_r._gm
            info.Lgm = self.profile_l._gm
            # info->Tgm was created when gm was still in local mode above
            # we cloned Tgm from the while profile was still locally configured, above
            status = libinfernal.cm_p7_modelconfig.p7_ProfileConfig5PrimeTrunc(info.Rgm, 100)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_ProfileConfig5PrimeTrunc")
            status = libinfernal.cm_p7_modelconfig.p7_ProfileConfig3PrimeTrunc(info.cm.fp7, info.Lgm, 100)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_ProfileConfig3PrimeTrunc")
            status = libinfernal.cm_p7_modelconfig.p7_ProfileConfig5PrimeAnd3PrimeTrunc(info.Tgm, 100)
            if status != libeasel.eslOK:
                raise UnexpectedError(status, "p7_ProfileConfig5PrimeAnd3PrimeTrunc")
        else:
            info.Rgm = NULL
            info.Lgm = NULL
            info.Tgm = NULL

        # copy E-value parameters
        libeasel.vec.esl_vec_FCopy(info.cm.fp7_evparam, libinfernal.CM_p7_NEVPARAM, info.p7_evparam)

        # compute msvdata
        info.msvdata = libhmmer.p7_scoredata.p7_hmm_ScoreDataCreate(info.om, NULL)
        if info.msvdata == NULL:
            raise AllocationError("P7_SCOREDATA", sizeof(P7_SCOREDATA))

        return 0

    cdef void free_info(
        self,
        WORKER_INFO *info,
    ) noexcept nogil:
        # TODO: free or use Python garbage collection with dedicated objects
        # CM_PIPELINE      *pli         # <-- owned by self
        # CM_TOPHITS       *th          # <-- owned by the returned Tophits
        # CM_t             *cm          # <-- owned by the input CM
        # P7_BG            *bg          # <-- owned by self.background
        # P7_OPROFILE      *om          # <-- owned by self.opt
        # P7_PROFILE       *gm          # <-- owned by self.profile
        # P7_PROFILE       *Rgm         # <-- owned by self.profile_r
        # P7_PROFILE       *Lgm         # <-- owned by self.profile_l
        # P7_PROFILE       *Tgm         # <-- owned by self.profile_t
        # P7_SCOREDATA     *msvdata     # <-- not owned
        # float            *p7_evparam  # <-- stack allocated (!)
        # float             smxsize
        if info.msvdata:
            libhmmer.p7_scoredata.p7_hmm_ScoreDataDestroy(info.msvdata)
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

        # ensure the CM defines a filter HMM
        if query.filter_hmm is None:
            raise ValueError(f"no filter HMM was found for CM {query.name!r}")

        # use struct to keep track of current worker state
        tinfo.p7_evparam = p7_evparam
        tinfo.smxsize = 128.0
        tinfo.cm = query._cm
        tinfo.pli = self._pli  # Maybe copy?
        tinfo.th = top_hits._th
        tinfo.bg = self.background._bg
        tinfo.Rgm = tinfo.Lgm = tinfo.Tgm = NULL
        tinfo.msvdata = NULL

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
        status = self._setup_hmm_filter(&tinfo, query)
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

        # Enforce threshold (and copy pipeline configuration) before returning
        top_hits._threshold(self)
        return top_hits


cdef class Alignment:
    cdef readonly Hit            hit
    cdef          CM_ALIDISPLAY* _ad

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Hit hit):
        self.hit = hit
        self._ad = hit._hit.ad

    def __str__(self):
        assert self._ad != NULL

        cdef int    status
        cdef object buffer = io.BytesIO()
        cdef FILE*  fp     = fopen_obj(buffer, "w")

        try:
            status = libinfernal.cm_alidisplay.cm_alidisplay_Print(
                fp,
                self._ad,
                0,
                -1,
                False,
            )
            if status == libeasel.eslEWRITE:
                raise OSError("Failed to write alignment")
            elif status != libeasel.eslOK:
                raise UnexpectedError(status, "cm_alidisplay_Print")
        finally:
            fclose(fp)

        return buffer.getvalue().decode("ascii")

    def __sizeof__(self):
        assert self._ad != NULL
        return sizeof(self) + libinfernal.cm_alidisplay.cm_alidisplay_Sizeof(self._ad)


    # --- Properties ---------------------------------------------------------

    @property
    def cm_from(self):
        """`int`: The start coordinate of the alignment in the query CM.
        """
        assert self._ad != NULL
        return self._ad.cfrom_emit

    @property
    def cm_to(self):
        """`int`: The end coordinate of the alignment in the query CM.
        """
        assert self._ad != NULL
        return self._ad.cto_emit

    @property
    def cm_name(self):
        """`str`: The name of the query CM.
        """
        assert self._ad != NULL
        assert self._ad.cmname != NULL
        return PyUnicode_FromString(self._ad.cmname)

    @property
    def cm_accession(self):
        """`str`: The accession of the query, or its name if it has none.
        """
        assert self._ad != NULL
        assert self._ad.cmacc != NULL
        return PyUnicode_FromString(self._ad.cmacc)

    @property
    def cm_description(self):
        """`str`: The description of the query, or its name if it has none.
        """
        assert self._ad != NULL
        assert self._ad.cmdesc != NULL
        return PyUnicode_FromString(self._ad.cmdesc)

    @property
    def cm_sequence(self):
        """`str`: The sequence of the query CM in the alignment.
        """
        assert self._ad != NULL
        return PyUnicode_DecodeASCII(self._ad.model, self._ad.N, NULL)

    @property
    def posterior_probabilities(self):
        """`str`: Posterior probability annotation of the alignment.
        """
        assert self._ad != NULL
        return PyUnicode_DecodeASCII(self._ad.ppline, self._ad.N, NULL)


    @property
    def target_from(self):
        """`int`: The start coordinate of the alignment in the target sequence.
        """
        assert self._ad != NULL
        return self._ad.sqfrom

    @property
    def target_name(self):
        """`str`: The name of the target sequence.

        .. versionchanged:: 0.12.0
            Property is now a `str` instead of `bytes`.

        """
        assert self._ad != NULL
        assert self._ad.sqname != NULL
        return PyUnicode_FromString(self._ad.sqname)

    @property
    def target_sequence(self):
        """`str`: The sequence of the target sequence in the alignment.
        """
        assert self._ad != NULL
        return PyUnicode_DecodeASCII(self._ad.aseq, self._ad.N, NULL)

    @property
    def target_to(self):
        """`int`: The end coordinate of the alignment in the target sequence.
        """
        assert self._ad != NULL
        return self._ad.sqto

    @property
    def identity_sequence(self):
        """`str`: The identity sequence between the query and the target.
        """
        assert self._ad != NULL
        return PyUnicode_DecodeASCII(self._ad.mline, self._ad.N, NULL)


cdef class Hit:
    # a reference to the TopHits that owns the wrapped CM_HIT, kept so that
    # the internal data is never deallocated before the Python class.
    cdef readonly TopHits   hits
    cdef readonly Alignment alignment
    cdef          CM_HIT*   _hit

    def __cinit__(self, TopHits hits, size_t index):
        assert hits._th != NULL
        assert index < hits._th.N
        self.hits = hits
        self._hit = hits._th.hit[index]
        self.alignment = Alignment(self)

    @property
    def name(self):
        """`str`: The name of the database hit.
        """
        assert self._hit != NULL
        assert self._hit.name != NULL
        return PyUnicode_FromString(self._hit.name)

    @name.setter
    def name(self, str name not None):
        assert self._hit != NULL

        cdef const char* data   = NULL
        cdef ssize_t     length = -1

        data = PyUnicode_AsUTF8AndSize(name, &length)
        self._hit.name = <char*> realloc(self._hit.name, max(1, sizeof(char) * length))
        if self._hit.name == NULL:
            raise AllocationError("char", sizeof(char), length)
        with nogil:
            memcpy(self._hit.name, data, sizeof(char) * (length + 1))

    @property
    def accession(self):
        """`str` or `None`: The accession of the database hit, if any.

        .. versionchanged:: 0.12.0
            Property is now a `str` instead of `bytes`.

        """
        assert self._hit != NULL
        if self._hit.acc == NULL:
            return None
        return PyUnicode_FromString(self._hit.acc)

    @accession.setter
    def accession(self, str accession):
        assert self._hit != NULL

        cdef const char* data   = NULL
        cdef ssize_t     length = -1

        if accession is None:
            free(self._hit.acc)
            self._hit.acc = NULL
        else:
            data = PyUnicode_AsUTF8AndSize(accession, &length)
            self._hit.acc = <char*> realloc(self._hit.acc, max(1, sizeof(char) * length))
            if self._hit.name == NULL:
                raise AllocationError("char", sizeof(char), length)
            with nogil:
                memcpy(self._hit.acc, data, sizeof(char) * (length + 1))

    @property
    def description(self):
        """`str` or `None`: The description of the database hit, if any.

        .. versionchanged:: 0.12.0
            Property is now a `str` instead of `bytes`.

        """
        assert self._hit != NULL
        if self._hit.desc == NULL:
            return None
        return PyUnicode_FromString(self._hit.desc)

    @description.setter
    def description(self, str description):
        assert self._hit != NULL

        cdef const char* data   = NULL
        cdef ssize_t     length = -1

        if description is None:
            free(self._hit.desc)
            self._hit.desc = NULL
        else:
            data = PyUnicode_AsUTF8AndSize(description, &length)
            self._hit.desc = <char*> realloc(self._hit.desc, max(1, sizeof(char) * length))
            if self._hit.name == NULL:
                raise AllocationError("char", sizeof(char), length)
            with nogil:
                memcpy(self._hit.desc, data, sizeof(char) * (length + 1))

    @property
    def score(self):
        """`float`: Bit score of the hit after correction.
        """
        assert self._hit != NULL
        return self._hit.score

    @property
    def bias(self):
        """`float`: The *null2*/*null3* contribution to the uncorrected score.
        """
        assert self._hit != NULL
        return self._hit.bias

    @property
    def pvalue(self):
        """`float`: The p-value of the bitscore.
        """
        assert self._hit != NULL
        return self._hit.pvalue

    @property
    def strand(self):
        """`str`: The strand where the hit is located (either "+" or "-").
        """
        assert self._hit != NULL
        return "+" if self._hit.start < self._hit.stop else "-"

    @property
    def included(self):
        """`bool`: Whether this hit is marked as *included*.
        """
        assert self._hit != NULL
        return self._hit.flags & libinfernal.cm_tophits.CM_HIT_IS_INCLUDED != 0

    @property
    def reported(self):
        """`bool`: Whether this hit is marked as *reported*.
        """
        assert self._hit != NULL
        return self._hit.flags & libinfernal.cm_tophits.CM_HIT_IS_REPORTED != 0


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

    @property
    def included(self):
        """iterator of `Hit`: An iterator over the hits marked as *included*.
        """
        return SizedIterator(
            self._th.nincluded,
            filter(operator.attrgetter("included"), self)
        )

    @property
    def reported(self):
        """iterator of `Hit`: An iterator over the hits marked as *reported*.
        """
        return SizedIterator(
            self._th.nreported,
            filter(operator.attrgetter("reported"), self)
        )

    # --- Utils --------------------------------------------------------------

    cdef int _threshold(self, Pipeline pipeline) except 1 nogil:
        # threshold the top hits with the given pipeline numbers
        cdef int status = libinfernal.cm_tophits.cm_tophits_Threshold(self._th, pipeline._pli)
        if status != libeasel.eslOK:
            raise UnexpectedError(status, "cm_tophits_Threshold")
        # record pipeline configuration
        # NOTE: the pointers on the copy are set to NULL by precaution, but since
        # `pipeline._pli` is allocated by the Python object directly, it will not
        # be deallocated, so there should be no risk of double free nevertheless.
        memcpy(&self._pli, pipeline._pli, sizeof(CM_PIPELINE))
        self._pli.oxf = self._pli.oxb = self._pli.fwd = self._pli.bck = NULL
        self._pli.gxf = self._pli.gxb = self._pli.gfwd = self._pli.gbck = NULL
        self._pli.r = NULL
        self._pli.ddef = NULL
        self._pli.cmfp = NULL
        return 0

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
        cdef str         sname  = None
        cdef str         sacc   = None
        cdef const char* qname  = NULL
        cdef const char* qacc   = NULL

        if isinstance(self._query, CM):
            qname = (<CM> self._query)._cm.name
            qacc  = (<CM> self._query)._cm.acc
        elif self._query is not None:
            if self._query.name is not None:
                sname = self._query.name
                qname = self._query.name
            if self._query.accession is not None:
                sacc = self._query.accession
                qacc = self._query.accession

        file = fopen_obj(fh, "w")
        try:
            if format == "3":
                fname = "cm_tophits_TabularTargets3"
                status = libinfernal.cm_tophits.cm_tophits_TabularTargets3(
                    file,
                    <char*> qname,
                    <char*> qacc,
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
            del sname
            del sacc


class NodeType(enum.IntEnum):
    #DUMMY = libinfernal.DUMMY_nd
    BIF  = libinfernal.BIF_nd
    MATP = libinfernal.MATP_nd
    MATL = libinfernal.MATL_nd
    MATR = libinfernal.MATR_nd
    BEGL = libinfernal.BEGL_nd
    BEGR = libinfernal.BEGR_nd
    ROOT = libinfernal.ROOT_nd
    END  = libinfernal.END_nd

class StateType:
    D  = libinfernal.D_st
    MP = libinfernal.MP_st
    ML = libinfernal.ML_st
    MR = libinfernal.MR_st
    IL = libinfernal.IL_st
    IR = libinfernal.IR_st
    S  = libinfernal.S_st
    E  = libinfernal.E_st
    B  = libinfernal.B_st
    EL = libinfernal.EL_st

# --- Module init code -------------------------------------------------------

init_ilogsum()
FLogsumInit()
p7_FLogsumInit()
