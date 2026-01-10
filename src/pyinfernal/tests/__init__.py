from . import (
    test_doctest,
    test_cm,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    suite.addTests(loader.loadTestsFromModule(test_cm))
    return suite
