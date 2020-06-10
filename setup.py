from setuptools import setup

setup(
    name="pispec",
    version="0.0",
    description="Prototype package for computing photon injection spectra",
    zip_safe=False,
    packages=["pispec"],
    author = 'Boris Bolliet',
    author_email = 'boris.bolliet@gmail.com',
    url = 'https://github.com/borisbolliet/pi_spec',
    download_url = 'https://github.com/borisbolliet/pi_spec/archive/master.zip',
    package_data={
        "pispec": ["ct_database/*"]
    },

)
