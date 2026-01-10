from libhmmer.p7_bg cimport P7_BG
from libhmmer.p7_profile cimport P7_PROFILE
from libhmmer.p7_scoredata cimport P7_SCOREDATA
from libinfernal.cm cimport CM_t
from libinfernal.cm_pipeline cimport CM_PIPELINE
from libinfernal.cm_tophits cimport CM_TOPHITS

if HMMER_IMPL == "VMX":
    from libhmmer.impl_vmx.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "SSE":
    from libhmmer.impl_sse.p7_oprofile cimport P7_OPROFILE
elif HMMER_IMPL == "NEON":
    from libhmmer.impl_neon.p7_oprofile cimport P7_OPROFILE

ctypedef struct WORKER_INFO:
    CM_PIPELINE      *pli
    CM_TOPHITS       *th
    CM_t             *cm
    P7_BG            *bg
    P7_OPROFILE      *om
    P7_PROFILE       *gm
    P7_PROFILE       *Rgm
    P7_PROFILE       *Lgm
    P7_PROFILE       *Tgm
    P7_SCOREDATA     *msvdata
    float            *p7_evparam
    float             smxsize