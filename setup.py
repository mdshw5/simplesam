from setuptools import setup
import os.path

setup(
    name='simplesam',
    use_scm_version={"local_scheme": "no-local-version"},
    setup_requires=['setuptools_scm'],
    author='Matthew Shirley',
    author_email='mdshw5@gmail.com',
    url='http://mattshirley.com',
    description='Simple pure Python SAM parser and objects for working with SAM records',
    license='MIT',
    classifiers=[
                "Development Status :: 5 - Production/Stable",
                "License :: OSI Approved :: MIT License",
                "Environment :: Console",
                "Intended Audience :: Science/Research",
                "Natural Language :: English",
                "Operating System :: Unix",
                "Programming Language :: Python :: 2.7",
                "Programming Language :: Python :: 3.10",
                "Topic :: Scientific/Engineering :: Bio-Informatics",
                ],
    py_modules=['simplesam'],
    scripts=['scripts/pileup.py'],
    install_requires=['six', 'ordereddict']
)