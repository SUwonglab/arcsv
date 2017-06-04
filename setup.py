from setuptools import setup, find_packages

with open('README.md') as f:
    arcsv_readme = f.read()

with open('LICENSE') as f:
    arcsv_license = f.read()

setup(
    name = 'ARC-SV',
    version = '0.9.0',
    description = 'Automated Reconstruction of Complex Structural Variants from WGS Data',
    long_description = arcsv_readme,
    author = 'Joey Arthur',
    author_email = 'josephgarthur@gmail.com',
    license = arcsv_license,
    # url = 'http://github.com/jgarthur/arcsv/',
    packages = find_packages(exclude = ('build','dist')),
    package_data = {'arcsv' : ['resources/*/*']},
    install_requires = [
        'python-igraph',
        'matplotlib',
        'numpy',
        'patsy',
        'pysam',
        'scikit-learn'
        # TODO need scipy?
    ],
    scripts = ['bin/arcsv']
    # TODO specific dependency links
)
