# mols2html
Pacote que permite gerar paginas html a partir de dataframes de moleculas, apresentar as suas propriedades, selecionar e descarregar em varios formatos, para acompanhar e apresentar resultados de variados modelos moleculares, sempre em grande estilo, com o visual ChEMBL.

# Instalação
```pip install --ignore-installed git+https://github.com/RaulSofia/mols2html.git```

# Utilização
Fornece um objeto [Mostrador()](/mols2html/mostrador.py) ao qual se podem ir adicionando DataFrames de dados moleculares bem como selecionando quais as colunas de propriedades que devem ser mostradas. Gera uma página HTML autocontida, sem links. Na pasta Um [exemplo](/examples/EXEMPLO.py) feito com este [dataset](/examples/dataset_teste.csv) gerou este [HTML](/examples/Resultados.html). Explicações mais aprofundadas em cada [função](/mols2html/mostrador.py).
