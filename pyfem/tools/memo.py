import functools

class memoize(object):
    def __init__(self, func):
        self.func = func
        self.memoized = {}
        self.method_cache = {}
    def __call__(self, *args):
        return self.cache_get(self.memoized, args,
            lambda: self.func(*args))
    def __get__(self, obj, objtype):
        #print id(self.method_cache)
        return self.cache_get(self.method_cache, obj,
            lambda: self.__class__(functools.partial(self.func, obj)))
    def cache_get(self, cache, key, func):
        try:
            return cache[key]
        except KeyError:
            cache[key] = func()
            return cache[key]
        except TypeError:
         # uncachable -- for instance, passing a list as an argument. 
         # Better to not cache than to blow up entirely. 
         return self.func(*args)

class Adder(object):
    @memoize
    def add(self, arg1, arg2):
        print 'CALCULATING', arg1, '+', arg2
        return arg1 + arg2


@memoize
def subtract(arg1, arg2):
    print 'CALCULATING', arg1, '-', arg2
    return arg1 - arg2

def main():
    print subtract(10, 5)
    print subtract(10, 5)

    adder1 = Adder()
    adder2 = Adder()
    print adder1.add(5, 5)
    print adder1.add(5, 5)
    print adder2.add(5, 5)
    print adder2.add(5, 3)

if __name__ == '__main__':
    main()
