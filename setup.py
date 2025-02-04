from setuptools import setup, find_packages

setup(
    packages=find_packages(),
    package_data={
        'polypal': ['data/*']
    },
    entry_points={
        'console_scripts': [
            'polypal=polypal.main:main',
        ]
    },
)