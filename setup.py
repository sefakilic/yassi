from distutils.core import setup, Extension

module1 = Extension('yassi',
                    sources = ['yassi.c'])

setup (name = 'yassi',
       version = '0.1',
       description = 'Yet Another Site Search Implementation',
       ext_modules = [module1])
