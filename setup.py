from setuptools import setup, find_packages


setup(
    name='ptes',
    version='0.1',
    packages=find_packages(exclude=find_packages("research").__add__(["ptes.tests", "research"])),
    package_data={'ptes': ['configs/*', 'tests/*']},
    install_requires=[
        'pandas>=0.23.4',
        'numpy>=1.15.1',
        'pyinterval>=1.2.0',
        'pyaml',
        'matplotlib>=2.2.2',
        ],
    entry_points={
            'console_scripts': [
                "ptes.random_shuffle = ptes.panhandles.random_shuffle:main",
                "ptes.bed_to_overlaps = ptes.panhandles.bed_to_overlaps:main",
                "ptes.star_chim_to_bed = ptes.panhandles.star_chim_to_bed:main",

            ]
      }
    )