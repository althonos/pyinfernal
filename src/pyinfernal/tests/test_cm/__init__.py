from . import (
    test_cmfile
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_cmfile))
    return suite
