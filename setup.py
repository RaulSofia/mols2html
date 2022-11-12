from setuptools import find_packages, setup

setup(
    name='mols2html',
    version='1.0.0',
    author='Raul Sofia',
    author_email='rauljcsofia.com',
    packages=['mols2html'],
    scripts=[],
    url='https://github.com/RaulSofia/mols2html',
    license='LICENSE.txt',
    description='Pacote que permite gerar paginas html a partir de dataframes de moleculas, apresentar as suas propriedades, selecionar e descarregar em varios formatos, para acompanhar e apresentar resultados de variados modelos moleculares, sempre em grande estilo, com o visual ChEMBL.',
    long_description=open('./README.txt').read(),
    # packages = find_packages(),
    install_requires=[
        "beautifulsoup4"
        "numpy"
        "pandas"
        "rdkit"
        "setuptools"
        "tqdm"
    ]
    
)
