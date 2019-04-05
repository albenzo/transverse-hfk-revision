from distutils.core import setup, Extension
import os

include = ['./src']
libs = []
lib_dirs = []

if os.name == 'posix':
    include += ['/usr/local/include']
    libs += ['argp']
    lib_dirs += ['/usr/local/lib']

_transverseHFK_module = Extension('tHFK._tHFK',
                                  sources = ['./tHFK/_transverseHFKmodule.c'],
                                  include_dirs = include,
                                  libraries = libs,
                                  library_dirs = lib_dirs)

setup(name = '_tHFK',
      version = '1.0',
      description = 'Computes transverse knot invariants',
      author = "Lucas Meyers",
      author_email = "lmeye22@lsu.edu",
      url = "https://github.com/albenzo/transverse-hfk-revision",
      packages = ['tHFK'],
      ext_modules = [_transverseHFK_module],
      long_description =
      """
      """)
