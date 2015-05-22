
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'version': '0.1',
    'description': 'Parsers',
    'author': 'mc',
    'install_requires': ['nose'],
    'packages': ['gtf'],
    'scripts': [],
    'name': 'gtf',
    'package_data': {'tests': ['data/*'],},
    'entry_points': {
        'console_scripts': [
            'parse-gtf = gtf.parse_gtf:main',
            ],
        },
}

setup(**config)
