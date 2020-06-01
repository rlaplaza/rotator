import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name="rotator",
    version="0.0",
    package_dir={"rotator": "rotator"},
    package=["rotator", "rotator/test"],
    author="R.LAPLAZA",
    author_email="laplazasolanas@gmail.com",
    description="Geometry manipulation of molecule files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rlaplaza/rotator",
    packages=setuptools.find_packages(),
    classifiers=["Programming Language :: Python :: 3"],
)
