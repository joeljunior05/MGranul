from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext
import os

from numpy.distutils.misc_util import get_numpy_include_dirs

COMPILER = "g++"

os.environ["CC"] = COMPILER
os.environ["CXX"] = COMPILER

include = get_numpy_include_dirs() + ['../cpp/lib', '../cpp/local/include']

msgranul = Extension('msgranul',
                    include_dirs    = include,
                    libraries = ['opencv_core',
                                    'opencv_imgproc',
                                    'opencv_features2d',
                                    'opencv_flann',
                                    'opencv_highgui'],
                    library_dirs= ['../cpp/local/lib'],
                    extra_compile_args  = ['-std=c++11', '-fno-stack-protector'],
                    extra_objects       = ['../cpp/local/lib/libgranul.a'],
                    sources             = ['msgranul.cpp']
                    )

class build_ext_subclass( build_ext ):
    def build_extensions(self):
        #c = self.compiler.compiler_type
        #print "compiler attr", self.compiler.__dict__
        #print "compiler", self.compiler.compiler
        #print "compiler is",c
        self.compiler.compiler([COMPILER])
        build_ext.build_extensions(self)



setup (name = 'MSGranul',
       version = '0.1',
       description = 'MSGranul wrapper',
       ext_modules = [msgranul])
