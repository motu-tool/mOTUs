from setuptools import setup
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

long_description = read('README.md')
setup(
    name = "motu-profiler",
    version = "3.0.3",
    author = "Alessio Milanese",
    author_email = "alessiom@ethz.ch",
    description = ("Taxonomic profiling of metagenomes from diverse environments with mOTUs3"),
    license = "GPLv3",
    include_package_data=True,
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords = "bioinformatics metagenomics taxonomic profiling",
    url = "https://github.com/motu-tool/mOTUs",
    packages=['motus'],
    download_url = "https://github.com/motu-tool/mOTUs/archive/refs/tags/3.0.3.tar.gz",
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    entry_points = {
        'console_scripts': ['motus=motus.motus:main'],
    },
)
