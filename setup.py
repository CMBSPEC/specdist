from setuptools import setup




setup(
    name="specdist",
    version="0.1",
    description="Python package for spectral distortions",
    zip_safe=False,
    packages=["specdist"],
    author = 'Boris Bolliet and Jens Chluba',
    author_email = 'boris.bolliet@gmail.com',
    url = 'https://github.com/borisbolliet/specdist',
    download_url = 'https://github.com/borisbolliet/specdist/archive/master.zip',
    package_data={
        "specdist": ["data/*txt"],
        #"data/ct_database/case_1_040520/*txt"]#,
    },


)
