import contextlib
import cProfile
import os


__all__ = ["profile", "profiler"]


def profile(path=None, quiet=False, *args, **kwargs):
    """
    Profile decorator.

    This function can be used as a decorator to profile the python
    code executed by the decorated function.

    .. code:: python

        import bio

        with bio.profiler():
            # Your code goes here

    Args:
        path (str, optional): File path for profiler output.
        quiet (bool): If true, print output. Otherwise do not.
        *args: Positional arguments passed to `cProfile.Profile`
        *args: Keyword arguments passed to `cProfile.Profile`
    """
    prof = cProfile.Profile(*args, **kwargs)

    def decorator(func):
        def wrapper(*fargs, **fkwargs):
            nonlocal prof, path, quiet
            prof.enable()
            result = func(*fargs, **fkwargs)
            prof.disable()

            if path:
                path = os.path.expandvars(os.path.expanduser(path))
                prof.dump_stats(path)

            if not quiet:
                prof.print_stats()

            return result

        return wrapper

    return decorator


@contextlib.contextmanager
def profiler(path=None, quiet=False, *args, **kwargs):
    """
    Profiler Context Manager.

    This function is a context manager which can be used to profile
    the python code within the scope of the context manager.

    .. code:: python

        import bio

        with bio.profiler():
            # Your code goes here

    Args:
        path (str): File path for profiler output.
        quiet (bool): If true, print output. Otherwise do not.
        *args: Positional arguments passed to `cProfile.Profile`
        *args: Keyword arguments passed to `cProfile.Profile`

    Returns:
        str: Path to profiler output file.
    """
    prof = cProfile.Profile(*args, **kwargs)
    prof.enable()
    yield
    prof.disable()

    if path:
        path = os.path.expandvars(os.path.expanduser(path))
        prof.dump_stats(path)

    if not quiet:
        prof.print_stats()
