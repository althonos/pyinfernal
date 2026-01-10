# coding: utf-8
# isort: skip_file
"""Cython bindings and Python interface to Infernal.
"""

# NOTE(@althonos): This needs to stay on top of every other import to ensure
#                  that PyHMMER dynamic libraries (`libeasel`, `libhmmer`) are
#                  loaded first and therefore added to the linker table: then
#                  when `pyinfernal.cm` is loaded, it will load `libinfernal`
#                  which links `libeasel` and `libhmmer`, but won't have to 
#                  resolve them because they will have been loaded already. 
from pyhmmer import errors, easel, plan7

from . import cm
from .cm import __version__

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__all__ = [
    "errors",
    "easel",
    "plan7",
    "cm",
]

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version of the
    library on `Read The Docs <https://pyhmmer.readthedocs.io/en/v{}/>`_.

    """.format(__version__)
