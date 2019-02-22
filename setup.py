from distutils.core import setup, Extension

_transverseHFK_module = Extension('transverseHFK._transverseHFK',
                                  sources = ['_transverse_module'],
                                  include_dirs = ['.'])

setup(name = '_transverseHFK',
      version = '1.0',
      description = 'Computes transverse knot invariants',
      packages = ['transverseHFK'],
      ext_modules = [_transverse_module])
