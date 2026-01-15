# coding: utf-8
# isort: skip_file
"""Cython bindings and Python interface to Infernal 1.1.

Infernal is a biological sequence analysis method that uses profile 
stochastic context-free grammars called *covariance models* (CMs) to 
identify RNA structure and sequence similarities. Infernal was developed 
by Eric P. Nawrocki during his PhD thesis in the 
`Eddy Laboratory <http://eddylab.org/>`_.

``pyinfernal`` is a Python package implemented with the
`Cython <https://cython.org/>`_ language that provides bindings to 
Infernal. It directly interacts with the Infernal internals. 
It builds on top of `pyhmmer <https://github.com/althonos/pyhmmer>`_
and follows a generally similar interface.

References:
    - Nawrocki, Eric P., Kolbe, Diana L., and Sean R. Eddy. 
      “Infernal 1.0: inference of RNA alignments.” 
      Bioinformatics (Oxford, England) vol. 25,10 (2009): 1335-7. 
      :doi:`10.1093/bioinformatics/btp157`. :pmid:`19307242`.
    - Nawrocki, Eric P., and Sean R. Eddy. 
      “Infernal 1.1: 100-fold faster RNA homology searches.” 
      Bioinformatics (Oxford, England) vol. 29,22 (2013): 2933-5. 
      :doi:`10.1093/bioinformatics/btt509`. :pmid:`24008419`.
    - Durbin, Richard, Eddy, Sean R., Krogh, Anders, and Graeme Mitchison.
      "RNA structure analysis". in *Biological sequence analysis: 
      Probabilistic models of proteins and nucleic acids*. 
      261–299 (Cambridge University Press, Cambridge, 1998). 
      :isbn:`978-0-511-79049-2`.

"""

# NOTE(@althonos): This needs to stay on top of every other import to ensure
#                  that PyHMMER dynamic libraries (`libeasel`, `libhmmer`) are
#                  loaded first and therefore added to the linker table: then
#                  when `pyinfernal.cm` is loaded, it will load `libinfernal`
#                  which links `libeasel` and `libhmmer`, but won't have to 
#                  resolve them because they will have been loaded already. 
from pyhmmer import errors, easel, plan7

from . import cm, infernal
from .cm import __version__
from .infernal import cmsearch

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"
__all__ = [
    "cm",
    "infernal",
    "cmsearch",
]

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version of the
    library on `Read The Docs <https://pyinfernal.readthedocs.io/en/v{}/>`_.

    """.format(__version__)
