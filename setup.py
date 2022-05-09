# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 17:22:56 2020

@author: amit
"""

from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("cython_node_val.pyx"),
    zip_safe=False,
)