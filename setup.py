#!/usr/bin/env python
# -*- encoding: utf8 -*-
#import glob
#import inspect
import io
import os

#from setuptools import find_packages
from setuptools import setup
from distutils.core import Extension
import numpy

long_description = """
Source code: https://github.com/aaspip/pyntfa""".strip() 


def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")).read()



ntfac_module = Extension('ntfacfun', sources=['pyntfa/src/main.c',
											  'pyntfa/src/ntfa_alloc.c',
											  'pyntfa/src/ntfa_blas.c',
											  'pyntfa/src/ntfa_divnnsc.c',
											  'pyntfa/src/ntfa_conjgrad.c',
											  'pyntfa/src/ntfa_weight2.c',
											  'pyntfa/src/ntfa_decart.c',
											  'pyntfa/src/ntfa_triangle.c',
											  'pyntfa/src/ntfa_ntriangle.c',
											  'pyntfa/src/ntfa_ntrianglen.c'	],
                                        include_dirs=[numpy.get_include()])

setup(
    name="pyntfa",
    version="0.0.2.0",
    license='MIT License',
    description="A python package of non-stationary time-frequency analysis for multi-dimensional multi-channel seismic data",
    long_description=long_description,
    author="pyntfa developing team",
    author_email="chenyk2016@gmail.com",
    url="https://github.com/aaspip/pyntfa",
    ext_modules=[ntfac_module],
    packages=['pyntfa'],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    keywords=[
        "seismology", "earthquake seismology", "exploration seismology", "array seismology", "denoising", "science", "engineering", "structure", "local slope", "filtering"
    ],
    install_requires=[
        "numpy", "scipy", "matplotlib", "pynpre"
    ],
    extras_require={
        "docs": ["sphinx", "ipython", "runipy"]
    }
)
