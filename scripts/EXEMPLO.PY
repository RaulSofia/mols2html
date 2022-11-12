from mols2html.mostrador import Mostrador
import pandas as pd


if __name__ == '__main__': #EXEMPLO
    disp = Mostrador() #cria um objeto
    data = pd.read_csv("./test_dataset.csv", sep=";") #le um ficheiro csv para pd.DataFrame
    disp.add(data, subtitle="AlogP") #armazena os conteudos da pd.DataFrame fornecida no objeto
    disp.add(data, subtitle=["AlogP", "Molecular Weight", "Aromatic Rings"]) #mesma coisa, mesmo dataset, para demonstrar que os dados podem ser repetidos, mas cada cartão mostrar campos diferentes
    # disp.add(data)
    # disp.add(data, title="AlogP", subtitle="Smiles")
    # print(disp.data[["_subtitulo", "_titulo"]])
    # print(data)
    # print(disp.data)
    disp.render()#sort_by="AlogP")
    disp.show()
    # disp.save("./resultados.html")
    # sleep(10)