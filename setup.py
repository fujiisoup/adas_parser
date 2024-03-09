#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
from setuptools import setup
import re

# load version form _version.py
VERSIONFILE = "adas_parser/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

# module

setup(
    name="adas_parser",
    version=verstr,
    author="Keisuke Fujii",
    author_email="fujiik@ornl.gov",
    description=("Python small library to download and read adas data"),
    license="BSD 3-clause",
    keywords="spectroscopy, atomic data",
    url="http://github.com/fujiisoup/adas_parser",
    include_package_data=True,
    ext_modules=[],
    packages=["adas_parser",],
    package_dir={"adas_parser": "adas_parser"},
    py_modules=["adas_parser.__init__"],
    test_suite="tests",
    install_requires="""
        numpy>=1.11
        scipy>=1.00
        """,
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
