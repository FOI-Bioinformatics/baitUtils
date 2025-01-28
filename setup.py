from setuptools import setup, find_packages

# Load version from baitUtils/_version.py
version = {}
with open("baitUtils/_version.py") as fp:
    exec(fp.read(), version)

setup(
    name='baitUtils',
    version=version['__version__'],  # Use dynamic version loading
    packages=find_packages(),
    install_requires=[
        'numpy>=1.21.0',
        'pandas>=1.3.0',
        'matplotlib>=3.4.0',
        'seaborn>=0.11.0',
        'scikit-learn>=0.24.0',
        'biopython>=1.78',
    ],
    entry_points={
        'console_scripts': [
            'baitUtils = baitUtils.__main__:main',
        ],
    },
    author='Andreas Sjödin',
    author_email='andreas.sjodin@gmail.com',
    description='baitUtils is a comprehensive toolkit for the analysis and visualization of bait sequences used in in-solution hybridization',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/FOI-Bioinformatics/baitUtils',  
    classifiers=[
        'Programming Language :: Python :: 3.12',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: MIT License',
    ],
    include_package_data=True,
    python_requires='>=3.11',
)