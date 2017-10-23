"""A utility module, containing various shared functions.
"""

import os

def get_script_dir():
    """Get the base folder of the script.
    
    Returns
    -------
    str
        The base folder of this script.
    """
    
    return os.path.dirname(os.path.realpath(__file__))

def get_dependency(*path_components):
    """Given a relative filename of a dependency, return the absolute path
    relative to the current script.
    """
    
    return os.path.join(*((get_script_dir(), "dependencies") + path_components))

def flatten(seq_of_seqs):
    """Flatten a sequence of sequences into a single sequence.
    
    Returns
    -------
    list
        The flattened list.
    """
    
    return [item for sublist in seq_of_seqs for item in sublist]

def make_dir(dir_path):
    """Makes sure the specified path exists as a folder.
    
    If the path doesn't exist, it recursively goes back in the folder tree
    and creates each subpath.
    
    Parameters
    ----------
    dir_path : str
        The folder to create (if necessary).
    """
    
    if not os.path.isdir(dir_path):
        make_dir(os.path.dirname(dir_path))
        os.mkdir(dir_path)
     
     
# This code was borrowed from https://wiki.python.org/moin/PythonDecoratorLibrary#Memoize

import collections
import functools

class memoized(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)