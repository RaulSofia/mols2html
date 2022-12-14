from setuptools import find_packages, setup

setup(
    name='mols2html',
    version='1.0.13',
    author='Raul Sofia',
    author_email='rauljcsofia.com',
    packages=["mols2html"],
    scripts=[],
    package_dir={"": "."},
    include_package_data=True,
    package_data={"mols2html": ["*.html"]},
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
