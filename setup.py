import os.path as op
from setuptools import setup


def readme(fname):
    """Utility function to read the README file."""
    return open(op.join(op.dirname(__file__), fname)).read()


class DownloadData(Command):
    """
    You can extend this function to download Dream Challenge and othe rdata at runtime.
    """
    user_options = []

    def initialize_options(self):
        """Abstract method that is required to be overwritten"""
        pass

    def run(self):
        """
        This is where the work happens.
        """
        pass

    def finalize_options(self):
        """Abstract method that is required to be overwritten"""
        pass


setup(
    name='az_dream',
    version='0.0.1'
    description='My submission for the AZ Dream Challene',
    author='strokach',
    author_email='alex.strokach@utoronto.ca',
    url='http://ostrokach.github.io',
    packages=['az_dream'],
    # package_data={'elaspic': ['data/*']},
    # long_description=read("README.rst"),
    # install_requires=meta['requirements']['run'],
    # tests_require=meta['test']['requires'],
    # entry_points={'console_scripts': meta['build']['entry_points']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Structural Biology",
        "Topic :: Bioinformatics",
    ],
    # cmdclass={'train': DownloadData},
)
