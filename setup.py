from setuptools import setup, find_packages
from os.path import abspath, dirname, join
import pathlib

# Fetches the content from README.md
# This will be used for the "long_description" field.
README_MD = open(join(dirname(abspath(__file__)), "README.md")).read()
PYTHON = ">=3.7"

here = pathlib.Path(__file__).parent

with open(here / "requirements.txt", "r") as f:
    REQUIRED = f.readlines()

setup(
    name="proppy",
    version="2.0.0",
    packages=find_packages(exclude=["test", "examples"]),

    # The description that will be shown on PyPI.
    # Keep it short and concise
    # This field is OPTIONAL
    description="PropPy is an open-source python software for propagating charged high-energy particles (cosmic rays) in astrophysical environments",

    # The content that will be shown on your project page.
    # In this case, we're displaying whatever is there in our README.md file
    # This field is OPTIONAL
    long_description=README_MD,
    long_description_content_type="text/markdown",
    url="https://gitlab.ruhr-uni-bochum.de/reichp2y/rwpropa",
    author="Patrick Reichherzer",
    author_email="patrick.reichherzer@ruhr-uni-bochum,de",
    

    # Classifiers help categorize your project.
    # For a complete list of classifiers, visit:
    # https://pypi.org/classifiers
    # This is OPTIONAL
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License"
    ],
    python_requires=PYTHON,
    install_requires=REQUIRED,

    # Keywords are tags that identify your project and help searching for it
    # This field is OPTIONAL
    keywords="astronomy, astrophysics",
)