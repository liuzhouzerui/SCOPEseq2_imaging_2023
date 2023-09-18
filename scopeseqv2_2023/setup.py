from setuptools import find_packages, setup

__version__ = '0.1.0'
requires = ['scipy',
            'numpy',
            'pandas',
            ]

setup(
    name='scopeseq',
    version=__version__,
    python_requires='>=3.6',
    install_requires=requires,
    author='Sims Lab',
    author_email='',
    description='SCOPE-seq V2, linking image Optical Barcode to sequencing Barcode',
    license="MIT",
    packages=find_packages(),
)
