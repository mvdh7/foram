import setuptools, foram

with open("README.md", "r") as fh:
    long_description = fh.read()
with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()
setuptools.setup(
    name="foram",
    version=foram.__version__,
    author=foram.__author__,
    author_email="m.p.humphreys@icloud.com",
    description="Kinetic model of foraminiferal calcification",
    url="https://github.com/mvdh7/foram",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
