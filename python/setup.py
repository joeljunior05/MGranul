from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext
import os

os.environ["CC"] = "clang"
os.environ["CXX"] = "clang"

msgranul = Extension('msgranul',
                    include_dirs    = ['../cpp/lib'],
                    libraries = ['opencv_core',
                                    'opencv_imgproc',
                                    'opencv_features2d',
                                    'opencv_stitching',
                                    'opencv_highgui'],
                    library_dirs= ['/usr/local/lib'],
                    extra_compile_args  = ['-std=c++11', '-fno-stack-protector'],
                    extra_objects       = ['../cpp/release/liblib_granul.a'],
                    sources             = ['msgranul.cpp']
                    )

class build_ext_subclass( build_ext ):
    def build_extensions(self):
        #c = self.compiler.compiler_type
        #print "compiler attr", self.compiler.__dict__
        #print "compiler", self.compiler.compiler
        #print "compiler is",c
        self.compiler.compiler(["clang"])
        build_ext.build_extensions(self)



setup (name = 'MSGranul',
       version = '0.1',
       description = 'MSGranul wrapper',
       ext_modules = [msgranul])