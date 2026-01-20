from . import (
    test_cmsearch
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_cmsearch))
    return suite
