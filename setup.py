from setuptools import setup, find_packages

setup(
    name='MAST',
    author='Samuel D. Lotz',
    author_email='samuel.lotz@salotz.info',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Numpy',
        'Pandas',
        'Scipy',
    ],
)
