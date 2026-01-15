import typing

# `typing.Literal`` is only available in Python 3.8 and later
if typing.TYPE_CHECKING:

    if sys.version_info >= (3, 8):
        from typing import Literal, TypedDict
    else:
        from typing_extensions import Literal, TypedDict  # type: ignore

    if sys.version_info >= (3, 11):
        from typing import Unpack
    else:
        from typing_extensions import Unpack

    from .plan7 import BIT_CUTOFFS, STRAND, Background
    from .easel import Alphabet

    BACKEND = Literal["threading", "multiprocessing"]
    PARALLEL = Literal["queries", "targets"]

    class PipelineOptions(TypedDict, total=False):
        alphabet: Alphabet
        background: typing.Optional[Background]
        seed: int
        Z: int
        E: float
        T: typing.Optional[float]
        incE: float
        incT: typing.Optional[float]
