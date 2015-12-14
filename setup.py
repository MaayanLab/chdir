"""Installs Python package of Characteristic Direction.
"""

from setuptools import setup, find_packages


setup(
    name='chdir',
    url='https://github.com/MaayanLab/chdir',
    author='Avi Ma\'ayan',
    author_email='avi.maayan@mssm.edu',
    version='1.0',
    packages=find_packages(),
    install_requires=['numpy']
)
