from distutils.core import setup, Extension

_transverseHFK_module = Extension('tHFK._tHFK',
                                  sources = ['./tHFK/_transverseHFKmodule.c'],
                                  include_dirs = ['./src'])

setup(name = '_tHFK',
      version = '1.0',
      description = 'Computes transverse knot invariants',
      packages = ['tHFK'],
      ext_modules = [_transverseHFK_module])
