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

def get_dependency(filename):
    """Given a relative filename of a dependency, return the absolute path
    relative to the current script.
    """
    
    return os.path.join(get_script_dir(), "dependencies", filename)

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