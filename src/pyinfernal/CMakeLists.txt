include(CythonExtension)

cython_extension(cm LINKS libinfernal DIRECTIVES -E HMMER_IMPL=${HMMER_IMPL})
