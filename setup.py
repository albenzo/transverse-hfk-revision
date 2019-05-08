from distutils.core import setup, Extension
import os

include = ['./src']
libs = []
lib_dirs = []

_transverseHFK_module = Extension('transHFK._transHFK',
                                  sources = ['./transHFK/_transverseHFKmodule.c', 'src/states.c', 'src/TransverseHFK.c'],
                                  include_dirs = include,
                                  libraries = libs,
                                  library_dirs = lib_dirs)

setup(name = '_transHFK',
      version = '1.0',
      description = 'Computes transverse knot invariants',
      author = "Lucas Meyers",
      author_email = "lmeye22@lsu.edu",
      url = "https://github.com/albenzo/transverse-hfk-revision",
      packages = ['transHFK'],
      ext_modules = [_transverseHFK_module],
      long_description =
      """
      """)
