# coding: utf-8
"""Reimplementation of Infernal binaries with the PyInfernal API.

Note:
    Functions of this module handle parallelization using threads to run
    searches in parallel for the different queries. If less queries are
    given, the number of threads will be reduced to avoid spawning idle
    threads.

"""

from ._cmsearch import cmsearch

__all__ = [
    "cmsearch",
]