import pytest


@pytest.fixture(scope="module", autouse=True)
def _benchmark_init_and_summarize(request):
    """
    For each test module under benchmarks/:
      - Call module.initialize() once before the first test runs.
      - Call module.summarize() once after all tests in that module complete.
    If either function is missing, do nothing.
    """
    mod = request.module
    init = getattr(mod, "initialize", None)
    if callable(init):
        init()

    # All tests in this module (including parametrized cases) run between yield.
    yield

    summ = getattr(mod, "summarize", None)
    if callable(summ):
        # To make sure printed output appears in the terminal,
        # run pytest with -s (disable output capture) or use the terminal reporter.
        summ()
