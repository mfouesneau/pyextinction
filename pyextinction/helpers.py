"""
This is a first collection of tools making the design easier
"""
import warnings
from .ezunits import unit, hasUnit


def isNestedInstance(obj, cl):
    """ Test for sub-classes types
    I could not find a universal test

    Parameters
    ----------
    obj: object instance
        object to test

    cl: Class
        top level class to test

    returns
    -------
    r: bool
        True if obj is indeed an instance or subclass instance of cl
    """
    tree = [ cl ]
    if hasattr(cl, '__subclasses'):
        for k in cl.__subclasses():
            if hasattr(k, '__subclasses'):
                tree += k.__subclasses__()
    return issubclass(obj.__class__, tuple(tree))


def type_checker(name, obj, tp):
    """ Check a given type and raise a type error if not correct

    Parameters
    ----------
    name: str
        name of the variable to show in the exception text

    obj: object
        object to check

    tp: type
        expected type of obj

    Raises
    ------
    :exc:TypeError:
        raises a TypeError if object is not of the correct type of a subclass of it
    """
    if not isNestedInstance(obj, tp):
        txt = 'Expected "{0:s}" of type {1:s}, got {2:s} instead.'
        raise TypeError(txt.format(name, str(tp.__name__), str(type(obj).__name__)))


def warning_on_one_line(message, category, filename, lineno, file=None,
                        line=None):
    return " {0:s}:{1:d} {2:s}:{3:s}".format(filename, lineno,
                                             category.__name__, str(message))


def missing_units_warning(name, defaultunit):
    """ Warn if any unit is missing

    Parameters
    ----------
    name: str
        name of the variable

    defaultunit: str
        default unit definition

    Raises
    ------
    warning: warnings.warn
        warn if units are assumed
    """
    warnings.formatwarning = warning_on_one_line
    msg = 'Variable {0:s} does not have explicit units. Assuming `{1:s}`\n'
    # stacklevel makes the correct code reference
    warnings.warn(msg.format(name, defaultunit), stacklevel=4)


def val_in_unit(varname, value, defaultunit):
    """ check units and convert to defaultunit or create the unit information

    Parameters
    ----------
    varname: str
        name of the variable

    value: value
        value of the variable, which may be unitless

    defaultunit: str
        default units is unitless

    Returns
    -------
    quantity: ezunits.Quantity
        value with units

    Example
    -------
    >>> r = 0.5
    >>> print(val_in_unit('r', r, 'degree'))
    # UserWarning: Variable r does not have explicit units. Assuming `degree`
    <Quantity(0.5, 'degree')>

    >>> r = 0.5 * unit['degree']
    >>> print(val_in_unit('r', r, 'degree'))
    <Quantity(0.5, 'degree')>
    """
    if not hasUnit(value):
        missing_units_warning(varname, defaultunit)
        return value * unit[defaultunit]
    else:
        return value.to(defaultunit)
