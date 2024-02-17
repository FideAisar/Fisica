"""Examples of decorators

See https://realpython.com/primer-on-python-decorators/

"""


def do_twice(func):
    """Run the decorated function twice"""

    def wrapper_do_twice(*args, **kwargs):
        func()
        func()

    return wrapper_do_twice
