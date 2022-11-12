from setuptools import find_packages, setup

setup(
    name='mols2html',
    version='1.0.1',
    author='Raul Sofia',
    author_email='rauljcsofia.com',
    packages=find_packages(where="mols2html"),
    scripts=[],
    package_dir={"": "mols2html"},
    include_package_data=True,
    url='https://github.com/RaulSofia/mols2html',
    license='LICENSE.txt',
    description='Pacote que permite gerar paginas html a partir de dataframes de moleculas, apresentar as suas propriedades, selecionar e descarregar em varios formatos, para acompanhar e apresentar resultados de variados modelos moleculares, sempre em grande estilo, com o visual ChEMBL.',
    long_description=open('README.md').read(),
    install_requires=[
        "beautifulsoup4",
        "numpy",
        "pandas",
        "rdkit",
        "setuptools",
        "tqdm",
    ]

)
