from setuptools import setup




setup(
    name="specdist",
    version="0.0",
    description="Prototype package for computing photon injection spectra",
    zip_safe=False,
    packages=["pispec"],
    author = 'Boris Bolliet and Jens Chluba',
    author_email = 'boris.bolliet@gmail.com',
    url = 'https://github.com/borisbolliet/pi_spec',
    download_url = 'https://github.com/borisbolliet/pi_spec/archive/master.zip',
    package_data={
        "specdist": ["data/*txt"],
        #"data/ct_database/case_1_040520/*txt"]#,
    },


)
