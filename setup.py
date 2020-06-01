import setuptools
import masp.version

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="masp",
    version=masp.version.get_API_version(),
    author="Andres Perez-Lopez",
    author_email="contact@andresperezlopez.com",
    url="https://github.com/andresperezlopez/masp",
    description="Multichannel Acoustic Signal Processing",
    long_description="MASP: a Multichannel Acoustic Signal Processing library for Python",
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)