PyInfernal |Stars|
==================

.. .. |Logo| image:: /_images/logo.png
..    :scale: 40%
..    :class: dark-light

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyinfernal.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyinfernal/stargazers
   :class: dark-light

*Cython bindings and Python interface to* `Infernal 1.1 <http://eddylab.org/infernal/>`_.

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|


.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyinfernal/test.yml?branch=master&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyinfernal/actions

.. |GitLabCI| image:: https://img.shields.io/gitlab/pipeline/larralde/pyinfernal/master?gitlab_url=https%3A%2F%2Fgit.embl.de&logo=gitlab&style=flat-square&maxAge=600
   :target: https://git.embl.de/larralde/pyinfernal/-/pipelines

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyinfernal?logo=codecov&style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyinfernal/

.. |PyPI| image:: https://img.shields.io/pypi/v/pyinfernal.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyinfernal

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyinfernal?ogo=anaconda&style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyinfernal

.. |AUR| image:: https://img.shields.io/aur/version/python-pyinfernal?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyinfernal

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyinfernal?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyinfernal/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyinfernal.svg?logo=python&style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyinfernal/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyinfernal.svg?logo=python&style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyinfernal/#files

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyinfernal/

.. |Mirror| image:: https://img.shields.io/badge/mirror-LUMC-003EAA.svg?maxAge=2678400&style=flat-square
   :target: https://git.lumc.nl/mflarralde/pyinfernal/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyinfernal.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyinfernal/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pyinfernal?style=flat-square&maxAge=3600
   :target: http://pyinfernal.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyinfernal/blob/master/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fpepy.tech%2Fprojects%2Fpyinfernal&search=%5B0-9%5D%2B.%5B0-9%5D%2B(k%7CM)&style=flat-square&label=downloads&color=303f9f&cacheSeconds=86400
   :target: https://pepy.tech/project/pyinfernal


Overview
--------

`Infernal <https://eddylab.org/infernal>`_ is a biological sequence analysis 
method that uses profile stochastic context-free grammars called 
*covariance models* (CMs) to identify RNA structure and sequence similarities. 
Infernal was developed by `Eric P. Nawrocki <https://scholar.google.com/citations?user=jn84-g0AAAAJ&hl=fr>`_
during his PhD thesis in the `Eddy Laboratory <http://eddylab.org/>`_.

``pyinfernal`` is a Python package, implemented using the `Cython <https://cython.org/>`_
language, that provides bindings to Infernal, directly interacting with the
Infernal internals. It builds on top of `PyHMMER <https://pyhmmer.readthedocs.io>`_, 
and follows a generally similar interface.

Setup
-----

Run ``pip install pyinfernal`` in a shell to download the latest release and all
its dependencies from PyPi, or have a look at the
:doc:`Installation page <guide/install>` to find other ways to install ``pyinfernal``.


Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   Examples <examples/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. include:: related.rst


License
-------

This library is provided under the `MIT License <https://choosealicense.com/licenses/mit/>`_.
The Infernal code is available under
the `BSD 3-clause <https://choosealicense.com/licenses/bsd-3-clause/>`_ license 
which allows redistribution of the sources in the ``pyinfernal`` 
distribution. See the :doc:`Copyright Notice <guide/copyright>` section 
for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `Infernal <http://eddylab.org/infernal>`_ *authors. It was developed by*
`Martin Larralde <https://github.com/althonos>`_ *during his PhD project
at the* `Leiden University Medical Center <https://www.lumc.nl/en/>`_
*in the* `Zeller team <https://zellerlab.org>`_.
