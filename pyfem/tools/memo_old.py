from numpy import ndarray
#def memoize(f):
#    cache= {}
#    def memf(*x):
#        y = list(x)
#        for i,arg in enumerate(y):
#            if isinstance(arg, ndarray) or isinstance(arg, list):
#                y[i] = tuple(arg)
#
#        y = tuple(y)
#        if y not in cache:
#            cache[y] = f(*x)
#        return cache[y]
#    return memf


def to_tuple(x):
    y = list(x)
    for i,arg in enumerate(y):
        if isinstance(arg, ndarray) or isinstance(arg, list):
            y[i] = tuple(arg)

    return tuple(y)

class memoize(object): 
    def __init__(self, func): 
        self.func = func 
        self.memoized = {} 
        self.method_cache = {} 
    def __call__(self, *args): 
        return self.cache_get(self.memoized, args, lambda: self.func(*args)) 
    def __get__(self, obj, objtype): 
        return self.cache_get(self.method_cache, obj, lambda: self.__class__(functools.partial(self.func, obj))) 
    def cache_get(self, cache, key, func): 
        key = to_tuple(key)
        try: 
            return cache[key] 
        except KeyError: 
            cache[key] = func() 
            return cache[key] 
        except TypeError: # unhashable argument
            return func() 

