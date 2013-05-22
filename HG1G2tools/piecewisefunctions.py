from numbers import Number

class UndefinedError(ArithmeticError):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        
class BoundsError(ArithmeticError):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class PiecewiseFunction:
    """Piecewise defined function."""
    def __init__(self):
        self.pieces = []
    def add_component(self, f, a, b):
        """Add a new component to piecewise function.
        
        f is a function, a and b are the lower and upper bounds."""
        # Check order of bounds
        if a>b:
            t = a
            a = b
            b = t
        elif a==b:
            raise BoundsError("Degenerate interval (%f,%f)." % (a,b))
        # Check for overlap with existing pieces
        for g,c,d in self.pieces:
            if (a<c and b>c) or (a>c and b<d) or (a<d and b>d):
                raise BoundsError("Overlapping intervals: (%f,%f) (%f,%f)" % (a,b,c,d))
        # Everything is ok, add to list
        self.pieces.append((f, a, b))
    def __call__(self, x):
        """Get value of piecewise function at x."""
        # Loop over pieces, they're guaranteed not to overlap.
        for f,a,b in self.pieces:
            if a<=x and x<=b:
                return f(x)
        raise UndefinedError("Piecewise function not defined at %f." % x)
    def __setitem__(self, key, value):
        key_is_tuple = (type(key)==tuple and len(key)==2)
        if not key_is_tuple:
            raise ValueError("Key must be a tuple of two numbers.")
        a = key[0]
        b = key[1]
        lower_not_number = not isinstance(a, Number)
        upper_not_number = not isinstance(b, Number)
        if lower_not_number or upper_not_number:
            raise ValueError("Key must be a tuple of two numbers.")
        # Check for correct type of value and try adding to pieces:
        if hasattr(value, '__call__'):
            self.add_component(value, a, b)
            return
        elif isinstance(value, Number):
            self.add_component(lambda x: value, a, b)
            return
        raise ValueError("Value must be a function or numeric constant.")
