from setuptools import setup, find_packages

with open('README.md') as f:
    arcsv_readme = f.read()

with open('LICENSE') as f:
    arcsv_license = f.read()

version = '0.0.0'
for line in open('arcsv/_version.py'):
    if line.startswith('__version__'):
        version = line.split('=')[1].strip().strip('"')

setup(
    name='arcsv',
    version=version,
    description='Automated Reconstruction of Complex Structural Variants from WGS Data',
    long_description=arcsv_readme,
    author='Joey Arthur',
    author_email='josephgarthur@gmail.com',
    license=arcsv_license,
    # url='http://github.com/jgarthur/arcsv/',
    packages=find_packages(exclude=('build', 'dist')),
    install_requires=[
        'pysam',
        'numpy',
        'scipy',
        'scikit-learn',
        'matplotlib',
        'patsy',
        'python-igraph',
        'pyinter'
    ],
    scripts=['bin/arcsv']
)
