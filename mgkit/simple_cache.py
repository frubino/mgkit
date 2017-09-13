class memoize(dict):
    """
    a cache found on the `PythonDecoratorLibrary <https://wiki.python.org/moin/PythonDecoratorLibrary#Alternate_memoize_as_dict_subclass>`_

    Not sure about the license for it.
    """
    def __init__(self, func):
        self.func = func

    def __call__(self, *args):
        return self[args]

    def __missing__(self, key):
        result = self[key] = self.func(*key)
        return result
