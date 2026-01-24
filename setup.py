from setuptools import setup, find_packages

setup(
    name="stix_flarelist_science",
    version="0.1.0",
    python_requires=">=3.11",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "stix-train=stix_train.cli:main",
            "stix-submit=stix_train.submit:main",
            "stix-fetch=generate_flarelist.cli:main",
            "stix-plot=stix_train.plot_results:main",
        ],
    },
    install_requires=[
        "stixpy @ git+https://github.com/TCDSolar/stixpy.git@main",
        "stixdcpy",
        "astrospice",
    ],
)
