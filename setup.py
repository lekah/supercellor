from setuptools import find_packages
from numpy.distutils.core import setup, Extension
from json import load as json_load

ext = Extension(
        name = 'supercellor.lib.optimal_supercell',
        sources = ['supercellor/lib/optimal_supercell.f90'],
    )
#~ ext2 = Extension(
        #~ name = 'samos.lib.mdutils',
        #~ sources = ['samos/lib/mdutils.f90'],
    #~ )

if __name__ == '__main__':
    with open('setup.json', 'r') as info:
        kwargs = json_load(info)
    setup(
        include_package_data=True,
        packages=find_packages(),
        package_data = {'': ['*.f90']},
        ext_modules = [ext],
        **kwargs
    )
