from setuptools import setup, find_packages

VERSION = '0.1.2' 
DESCRIPTION = 'Package for solving quantum heat problems in the linear regime'
LONG_DESCRIPTION = ''
print(find_packages())
# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="picoQuantum", 
        version=VERSION,
        author="Andrew Guthrie",
        author_email="<andrew.guthrie@aalto.fi>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'first package'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)


