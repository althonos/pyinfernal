# 🐍 PyInfernal [![Stars](https://img.shields.io/github/stars/althonos/pyinfernal.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyinfernal/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [Infernal](http://eddylab.org/infernal).*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyinfernal/test.yml?branch=master&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyinfernal/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyinfernal?logo=codecov&style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/pyinfernal/)
[![PyPI](https://img.shields.io/pypi/v/pyinfernal.svg?logo=pypi&style=flat-square&maxAge=3600)](https://pypi.org/project/pyinfernal)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyinfernal?logo=anaconda&style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/pyhmmer)
[![AUR](https://img.shields.io/aur/version/python-pyinfernal?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyinfernal)
[![Wheel](https://img.shields.io/pypi/wheel/pyinfernal.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyinfernal/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyinfernal.svg?logo=python&style=flat-square&maxAge=3600)](https://pypi.org/project/pyinfernal/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyinfernal.svg?logo=python&style=flat-square&maxAge=3600&label=impl)](https://pypi.org/project/pyinfernal/#files)
[![License](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/mit/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyinfernal/)
[![Mirror](https://img.shields.io/badge/mirror-LUMC-003EAA.svg?maxAge=2678400&style=flat-square)](https://git.lumc.nl/mflarralde/pyinfernal/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyinfernal.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyinfernal/issues)
[![Docs](https://img.shields.io/readthedocs/pyinfernal/latest?style=flat-square&maxAge=600)](https://pyinfernal.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyinfernal/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fpepy.tech%2Fprojects%2Fpyinfernal&search=%5B0-9%5D%2B.%5B0-9%5D%2B(k%7CM)&style=flat-square&label=downloads&color=303f9f&cacheSeconds=86400)](https://pepy.tech/project/pyinfernal)


## 🗺️ Overview

[Infernal](https://eddylab.org/infernal) is a biological sequence analysis method that uses profile stochastic 
context-free grammars called *covariance models* (CMs) to identify RNA structure and sequence similarities. 
Infernal was developed by [Eric P. Nawrocki](https://scholar.google.com/citations?user=jn84-g0AAAAJ&hl=fr) 
during his PhD thesis in the [Eddy Laboratory](http://eddylab.org/).

`pyinfernal` is a Python package, implemented using the [Cython](https://cython.org/)
language, that provides bindings to Infernal. It directly interacts with the
Infernal internals. It builds on top of [`pyhmmer`](https://github.com/althonos/pyhmmer)
and follows a generally similar interface.

*This library is still very experimental and has not been thoroughly tested yet, use with caution*.

## 🔧 Installing

`pyinfernal` can be installed from [PyPI](https://pypi.org/project/pyinfernal/),
which hosts some pre-built CPython wheels for Linux and MacOS on x86-64 and Arm64, as well as the code required to compile from source with Cython:
```console
$ pip install pyinfernal
```

Compilation for UNIX PowerPC is not tested in CI, but should work out of the
box. Note than non-UNIX operating systems (such as Windows) are not
supported by Infernal.

## 🔖 Citation

PyInfernal is scientific software, and builds on top of Infernal. Please cite
the [Infernal 1.1 application note](https://academic.oup.com/bioinformatics/article/29/22/2933/316439) 
in Bioinformatics, for instance:

> PyInfernal, a Python library binding to Infernal (Nawrocki & Eddy, 2013).

<!-- Detailed references are available on the [Publications page](https://pyinfernal.readthedocs.io/en/stable/guide/publications.html) of the
[online documentation](https://pyinfernal.readthedocs.io/). -->

Also refer to the [Infernal User's Guide](http://eddylab.org/infernal/Userguide.pdf)
which contains a section about citation and reproducibility.

<!-- ## 📖 Documentation

A complete [API reference](https://pyinfernal.readthedocs.io/en/stable/api/) can
be found in the [online documentation](https://pyinfernal.readthedocs.io/), or
directly from the command line using
[`pydoc`](https://docs.python.org/3/library/pydoc.html):
```console
$ pydoc pyinfernal.cm
$ pydoc pyinfernal.infernal
``` -->

## 💡 Example

Use `pyinfernal` to run `cmsearch` to search for the genome of 
*Escherichia coli str. K-12 substr. MG1655* ([U00096.3](www.ncbi.nlm.nih.gov/nuccore/U00096.3))
for models from [RFam](https://rfam.org). This will produce an iterable 
over [`TopHits`] that can be used for further sorting/querying in Python. 
Processing happens in parallel using Python threads, 
and a [`TopHits`] object is yielded for every [`CM`] in the input iterable.

[`CM`]: https://pyinfernal.readthedocs.io/en/stable/api/cm/hmms.html#pyinfernal.plan7.HMM
[`TopHits`]: https://pyinfernal.readthedocs.io/en/stable/api/cm/results.html#pyinfernal.plan7.TopHits

```python
import pyhmmer.easel
import pyinfernal.cm
import pyinfernal.infernal

with pyhmmer.easel.SequenceFile("U00096.3.faa", digital=True) as seq_file:
    sequences = seq_file.read_block()

with pyinfernal.cm.CMFile("RFam.cm") as cm_file:
    for hits in pyinfernal.cmsearch(cm_file, sequences, cpus=4):
        print(f"CM {hits.query.name.decode()} found {len(hits)} hits in the target sequences")
```

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/pyinfernal/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See [`CONTRIBUTING.md`](https://github.com/althonos/pyinfernal/blob/master/CONTRIBUTING.md) for more details.


<!-- ## 🔍 See Also -->


## ⚖️ License

This library is provided under the [MIT License](https://choosealicense.com/licenses/mit/).
The Infernal code is available under the
[BSD 3-clause](https://choosealicense.com/licenses/bsd-3-clause/) license.
See `vendor/infernal/LICENSE` for more information.

*This project is in no way affiliated, sponsored, or otherwise endorsed by
the [original Infernal authors](http://eddylab.org/infernal/). It was developed by
[Martin Larralde](https://github.com/althonos/pyinfernal) during his PhD project
at the [Leiden University Medical Center](https://www.lumc.nl/en/) in
the [Zeller Lab](https://zellerlab.org/).*
