from setuptools import setup, find_packages
from os.path import abspath, dirname, join

# Fetches the content from README.md
# This will be used for the "long_description" field.
README_MD = open(join(dirname(abspath(__file__)), "README.md")).read()

setup(
    name="proppy",
    version="2.0.0",
    packages=find_packages(exclude=["test", "examples"]),

    # The description that will be shown on PyPI.
    # Keep it short and concise
    # This field is OPTIONAL
    description="PropPy is an open-source python software for propagating charged high-energy particles (cosmic rays) in a astrophysical environments",

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
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Astronomers",
    ],

    # Keywords are tags that identify your project and help searching for it
    # This field is OPTIONAL
    keywords="astronomy, astrophysics",
)