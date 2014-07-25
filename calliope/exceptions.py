"""
Copyright (C) 2013 Stefan Pfenninger.
Licensed under the Apache 2.0 License (see LICENSE file).

exceptions.py
~~~~~~~~~~~~~

Exceptions and Warnings.

"""


class ModelError(Exception):
    """
    ModelErrors should stop execution of the model, e.g. due to a problem
    with the model formulation or input data.

    """
    pass


class OptionNotSetError(ModelError):
    pass


class ModelWarning(Warning):
    """
    ModelWarnings should be raised for errors where execution could still
    continue, i.e. they may be caught and suppressed in some contexts.

    """
    pass