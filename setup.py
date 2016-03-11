# based on https://github.com/pypa/sampleproject
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='mergevcf',

    version='1.0.1',

    description='Merge VCF calls',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/ljdursi/mergevcf',

    # Author details
    author='Jonathan Dursi',
    author_email='Jonathan.Dursi@oicr.on.ca',

    # Choose your license
    license='GPL',

    classifiers=[
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 2.8',
#        'Programming Language :: Python :: 3',
#        'Programming Language :: Python :: 3.2',
#        'Programming Language :: Python :: 3.3',
#        'Programming Language :: Python :: 3.4',
    ],

    keywords='merge vcfs',

    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    install_requires=['pyvcf'],

    test_suite='tests',

    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage'],
    },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
#    package_data={
#        'sample': ['package_data.dat'],
#    },

    entry_points={
        'console_scripts': [
            'mergevcf=mergevcf:main',
        ],
    },
)
