from setuptools import setup
import sys,os

with open('snipgenie/description.txt') as f:
    long_description = f.read()

setup(
    name = 'snipgenie',
    version = '0.6.0',
    description = 'variant calling and phylogenies from microbial WGS data',
    long_description = long_description,
    url='https://github.com/dmnfarrell/snipgenie',
    license='GPL v3',
    author = 'Damien Farrell',
    author_email = 'farrell.damien@gmail.com',
    packages = ['snipgenie'],
    package_data={'snipgenie': ['data/*.*','logo.png','description.txt',
                                  'styles/*.qss','icons/*.png',
                                  'plugins/*.py','plugins/icons/*.png',
                  'description.txt']
                 },
    install_requires=['numpy>=1.2',
                      'pandas>=0.24',
                      'matplotlib>=3.0',
                      'biopython>=1.5',
                      'pyvcf3',
                      'pyfaidx',
                      'bcbio_gff',
                      #'pyside2>=5.1',
                      #'toytree',
                      ],
    entry_points = {
        'console_scripts': [
            'snipgenie-gui=snipgenie.gui:main',
            'snipgenie=snipgenie.app:main',
            'snipgenie-treeview=snipgenie.treeview:main']
            },
    classifiers = ['Operating System :: OS Independent',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.7',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords = ['bioinformatics','biology','genomics']
)
