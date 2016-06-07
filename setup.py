from setuptools import setup
import os.path


def get_version(string):
    """ Parse the version number variable __version__ from a script. """
    import re
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str

setup(
    name='simplesam',
    version=get_version(open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'simplesam.py')).read()),
    author='Matthew Shirley',
    author_email='mdshw5@gmail.com',
    url='http://mattshirley.com',
    description='Simple pure Python SAM parser and objects for working with SAM records',
    license='MIT',
    classifiers=[
                "Development Status :: 3 - Alpha",
                "License :: OSI Approved :: MIT License",
                "Environment :: Console",
                "Intended Audience :: Science/Research",
                "Natural Language :: English",
                "Operating System :: Unix",
                "Programming Language :: Python :: 2.7",
                "Programming Language :: Python :: 3.4",
                "Topic :: Scientific/Engineering :: Bio-Informatics",
                ],
    py_modules=['simplesam'],
    install_requires=['six']
)
