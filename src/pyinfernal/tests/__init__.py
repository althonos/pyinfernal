from . import (
    test_doctest,
    test_cm,
    test_infernal,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    suite.addTests(loader.loadTestsFromModule(test_cm))
    suite.addTests(loader.loadTestsFromModule(test_infernal))
    return suite
