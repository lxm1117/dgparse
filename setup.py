from setuptools import setup, find_packages


setup(
    name='genomebrowser',
    version='0.0.3d',
    packages=find_packages(),
    install_requires=[
        'pytest',
    ],
    entry_points={
        "distutils.commands": ["run = genomebrowser:run"],
        'console_scripts': ['genomebrowser = genomebrowser:run'],
        'paste.app_factory': ['main = genomebrowser:main'],
    },
)
