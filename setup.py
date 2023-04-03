from setuptools import setup, find_packages

VERSION = 'v1.5.0'
DESCRIPTION = 'Electric field analyses for molecular simulations.'

# Setting up
setup(
    name="tupa",
    version=VERSION,
    author="Marcelo D. Poleto, Justin A. Lemkul",
    author_email="Marcelo D. Poleto <mdpoleto@vt.edu>",
    description=DESCRIPTION,
    packages=find_packages(),

    install_requires=['MDAnalysis>=2.2.0','numpy>1.2', 'configparser', 'argparse'],
    extras_require = {},
    keywords=['python', 'molecular dynamics', 'electric fields', 'MD trajectory analysis'],
    url="https://github.com/mdpoleto/tupa/",
    classifiers=[
        "Programming Language :: Python :: 3.x",
        "License :: OSI Approved :: GNU General Public License v3",
    ],
)
