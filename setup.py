from setuptools import setup

setup(
    name="pispec",
    version="0.0",
    description="Prototype package for computing photon injection spectra",
    zip_safe=False,
    packages=["pispec"],
    package_data={
        "pispec": ["ct_database/*"]
    },

)
