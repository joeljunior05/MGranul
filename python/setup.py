from distutils.core import setup, Extension

msgranul = Extension('msgranul',
                    include_dirs    = ['/usr/local/lib/python3.6/site-packages/numpy/core/include',
                                        '/usr/local/Cellar/opencv@2/2.4.13.6_2/include', 
                                        '../cpp/lib'],
                    libraries = ['opencv_core',
                                    'opencv_imgproc',
                                    'opencv_features2d'],
                    library_dirs= ['/usr/local/Cellar/opencv@2/2.4.13.6_2/lib'],
                    extra_compile_args  = ['-std=c++11'],
                    extra_objects       = ['../cpp/release/liblib_granul.a'],
                    sources             = ['msgranul.cpp']
                    )

setup (name = 'MSGranul',
       version = '0.1',
       description = 'MSGranul wrapper',
       ext_modules = [msgranul])