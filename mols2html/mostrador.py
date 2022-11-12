import base64
import os
import shutil
import time
import webbrowser
from copy import copy
from io import BytesIO
from time import sleep
from tkinter.filedialog import asksaveasfilename
import numpy as np

import pandas as pd
from bs4 import BeautifulSoup as bs

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdCoordGen import AddCoords
from tqdm import tqdm

package_dir = os.path.dirname(os.path.abspath(__file__))

class Mostrador():
    def __init__(self) -> None:
        """Construtor"""
        self.data = None
        self.tempdir = os.path.join(package_dir, "temp", str(int(time.time())))
        os.makedirs(self.tempdir)


    def add(self, data: pd.DataFrame, smiles="Smiles", title="Smiles", subtitle=None):
        """Adiciona dados a serem representados em html a uma instância de Mostrador()"""
        aux = data.copy(deep=True)
        # checkar se tem as colunas necessárias
        if not smiles or smiles not in data:
            print("Erro a adicionar ao mostrador: Falta a coluna", smiles)
            return
        if title:
            if title not in data:
                print("Erro a adicionar ao mostrador: Falta a coluna", title)
                return
            aux["_titulo"] = data[title].map(str)
        else:
            aux["_titulo"] = np.nan
        if subtitle:
             #usar este sep para ter a certeza q nao ha conflitos
            if type(subtitle) == str:
                subtitle = [subtitle]
            
            aux["_subtitulo"] = subtitle[0] + ":--->:" + data[subtitle[0]].map(str)
            for sub in subtitle[1:]:
                if sub not in data:
                    print("Erro a adicionar ao mostrador: Falta a coluna", sub)
                    return
                aux["_subtitulo"] += "|||" + sub + ":--->:" + data[sub].map(str)
        else:
            aux["_subtitulo"] = np.nan
            

        aux["Smiles"] = data[smiles] #normalizar os datasets
        aux["_hash"] = data[smiles].apply(hash) #------ ATENÇÃO ------- a func hash() nao mantem os mesmos valores entre execuções de python, é randomizada cada vez que o interpretador é inicializado, mas como aqui so esta a aser usada para substituir os carateres estranhos dos smiles por nums nao ha problema
        if self.data is None:
            self.data = aux
        else:
            self.data = pd.concat([self.data, aux], ignore_index=True)


    def clean(self):
        """Limpa dados do objeto em causa"""
        self.data = None
        index_file = os.path.join(self.tempdir, "index.html")
        if os.path.exists(index_file):
            os.remove(index_file)

    def render(self, sort_by=None, name="Resultados"):
        """Computa a página HTML, preenche e por aí fora, mas não a apresenta, fica armazenada na pasta temp, até à morte do objeto. Aceita um argumento sort_by que identifica a coluna pela qual o ordenamento deve ser feito, ou None por default"""
        if sort_by:
            self.data = self.data.sort_values(sort_by)
        print("A renderizar index.html...")

        #preencher html
        with open(os.path.join(package_dir, "index.html"), "r") as template_file:
            template = bs(template_file, "html.parser")
        
        template.title.string = name
        
        with open(os.path.join(package_dir, "cartao_molecula_elemento.html"), "r") as elemento_molecula:
            e_molecula = bs(elemento_molecula, "html.parser")

        div_moleculas = template.find(id="div-moleculas")
        print("A popular index.html...")
        template.find(id = "n_compostos").string = str(self.data.shape[0]) + " Compounds"

        #COISAS Q SE FAZEM P CADA MOLECULA
        for i, row in tqdm(self.data.iterrows(), total=self.data.shape[0]):
            cartao_molecula = copy(e_molecula)
            main_div = cartao_molecula.find()
            main_div["id"] = i
            main_div["data-details"] = row.to_json()
            #cartao_molecula.find(id="img-molecula")['src'] = os.path.join("imagens", str(row["_hash"]) + ".png") # ATENÇÃOOOOO, QUANDO O INDEX.HTML CORRE NO BROWSER JA NAO ESTAS COM O MESMO CWD, PELO Q O ./ REMETE PARA A PASTA TEMP/QQCOISA---------------------------------------------OLHAAA
            encoded_string = Mostrador.__render_image(row["Smiles"])
            if encoded_string:
                cartao_molecula.find(id="img-molecula")['src'] = "data:image/png;base64," + encoded_string
                cartao_molecula.find(id="label-id-molecula").string = str(i + 1)
                tag0, tag1 = cartao_molecula.find_all(id="BCK-cards-ID_MOLECULA-select")
                tag0["id"] = "BCK-cards-" + str(i) + "-select"
                tag1["for"] = "BCK-cards-" + str(i) + "-select"

            if pd.notna(row["_titulo"]):
                cartao_molecula.find(id="titulo").string = row["_titulo"]
            if pd.notna(row["_subtitulo"]):
                parent_subtitulo_elem = cartao_molecula.find(id="subtitulos")
                for s in row["_subtitulo"].split("|||"):
                    subtitulo_elem = cartao_molecula.new_tag("p")
                    t_subt, str_subt = s.split(":--->:")
                    t_subt_tag = cartao_molecula.new_tag("b")
                    t_subt_tag.string = t_subt + ": "
                    subtitulo_elem["class"] = "p-oneline"
                    subtitulo_elem.string = str_subt
                    subtitulo_elem.insert(0, t_subt_tag)
                    parent_subtitulo_elem.append(subtitulo_elem)
            div_moleculas.append(cartao_molecula)
    
        with open(os.path.join(self.tempdir, "index.html"), "w") as out:
            out.write(str(template))



    def show(self, sort_by=None, name="Resultados"):
        """Apresenta a página anteriormente renderizada. Caso não tenha ainda sido renderizada, chama self.render() para o fazer, e depois abre no browser"""
        index_file = os.path.join(self.tempdir, "index.html")
        if not os.path.exists(index_file): #redundancia para o caso de me esquecer
            self.render(sort_by, name)
        webbrowser.open(index_file)
        sleep(2) #como isto abre threads q n controlo o melhor é usar um sleep para dar tempo

    def __render_image(smile):
        """Calcula os gráficos 2d para o smiles fornecido"""
        try:
            mol = Chem.MolFromSmiles(smile)
        except:
            return "iVBORw0KGgoAAAANSUhEUgAAA0gAAAI8CAIAAABEQRzvAAAACXBIWXMAAA7EAAAOxAGVKw4bAAAgAElEQVR4nO3deZgU1dn+8XOqeplphoHGARVRQRCJYAwqKCqggkajRlHzQtBo3DCuLy7RGJf4atxXRFyjaIwLP0FRQCBscUkEVBQVowKCiAgitCwzTE931fn9MSRuzNbbU8v3c+XyQqT6ua8w3X33qa5T2hijcpVKpXI+NplMMpe5zGUuc5nLXOYyt4BzrZyPBAAAgKdQ7AAAAAKCYgcAABAQFDsAAICAoNgBAAAEBMUOAAAgICh2AAAAAUGxAwAACAiKHQAAQEBQ7AAAAAKCYgcAABAQFDsAAICAoNgBAAAEBMUOAAAgICh2AAAAAUGxAwAACAiKHQAAQEBQ7AAAAAKCYgcAABAQkVQqlfPByWQy52OZy1zmMpe5zGUuc5lb2Lms2AEAAAQExQ4AACAgKHYAAAABQbEDAAAICIodAABAQFDsAAAAAoJiBwAAEBAUOwAAgICg2AEAAARERDoAADQlmzWplFm3zv7iC715s9q8WW/erKur9ebNqqZG19Soujqrtlal0zqdNpmMrq1VrqscRzlO/eGWUlmljFIqElFaK9ve+r94XMfjqqxM/fefiYSuqFCVlaqiQrVurSsqVOvWlm2btm1N27bKtoX/rwCARlHsAHiC+eYbs3q1+uors2aNWrPGrFmj1q5V33xj1q/XmzbV/5lEfiO0UiqbVUqpTKaxJD/6dav//KvbqpVJJt1k0iSTpkMHU1Xltm9f/wvTpk1+6QCgACh2AErLdd3Vq9XKlWblyq3//Pxz88UXuoGypUscr1FWdbWqrrZXrvzxf3IjEdOxo7vjjqZjR6djx62/7tCh9CEBhBnFDkAxua5etcpesUKvWGEtX26vWJH5/HNdf4b0OzzV3nJjZbNqxQp7xYrv/qYbiWQ7d1ZduujddtP1/9xpJ2Xx5WYAxUKxA1BQdXXWsmXW0qX2kiX2kiV6+XKr/uxnKFnZrFqyRC1ZYv5zVtdEo3r33XX37qp7d92jh7XbbioeF04JIEAodgDy47p6+XL744/tjz6yP/5Yr1jBelQjdCajPvzQfPihUsoo5Silu3TRPXuqXr10r15Wly6s5wHIB8UOQMtVV9uLFtmLFlkffWR/8omVTksH8iutlFq2zCxbpiZPNko5ZWW6Z0+9555677313ntLpwPgPxQ7AM2TSrnvvWfeece8/XarpUstY5o+BC2ka2vV22+bt982Tz5ptE7stpvTq5ez115Or16mslI6HQAfoNgBaNiWLe5775l588ybb6olS/7725wsLAFtjL10qb10qXrxRaWU06WL07t3dt99nZ49VSwmnQ6AR1HsAPyQu3ixmTvXzJ1rFi788RWsEGEvW2YvWxZ7/nk3EnF69nT22Sfbp4/p3Fk6FwBv0SaP8ympVCrnY5PJJHOZy1wPza2rs999NzJ/vj1/vv311zk/LErJdOhg9e+vDzrI2mefH1xd65WfK+Yyl7mlncuKHRBqesMGe+7cyNy59rvvcg2E7+ivvjITJpgJE5x4XPftqwcMsPr319wDAwgxih0QRubrr6Mvv2z/61/2woV8YS4AdDqtXnvNvPZaVim933760EN1796mbVvpXABKjWIHhIhZv96dNcvMmmUWLiyTDoNi0Eqpt94yb73VSilnzz2zAwdm+/en4QHhQbEDgs9s3mxefdWdPt28+aY2RgXiFl5onFYq8uGHkQ8/dB94wNlnH+fQQzMHHqjKy6VzASguih0QXNmsO3eu+/LL5vXXdSaj6HOhZCllLVgQXbAgNmpU9oADsoMHO/vuq2xbOheAoqDYAQHkfvqpmTLFnTpVp1KKPgellFJWNht7/fXY6687yWR20KDMz39udtpJOhSAAqPYAcFhNm92Z840kyapDz9U9Dk0wE6l7PHj4+PHO3vsUXfEEdlDDuEULRAYFDsgCNzFi82ECe706bq2VjoLfMP++OPyjz92H344O2hQ3THHsN0xEAAUO8DP0ml39mx3wgS1aJFiiQ45sdLp2Msvx15+2enRo+7YY7MHHSSdCEDuKHaAL+l166KTJmWmTdMbN0pnQUDYH31U/tFH7kMPOSeeaA0Zotu3l04EoMUodoDPWB9/HH3hhchrr1l53A8QaIi1caMZOzb717/qww6zhg2z9txTOhGAFqDYAT7hupF//Ss2YYL98cfSURB82nHUjBnujBlur17W8OHWwIHK4h4lgA9Q7ADPq6uLzJwZHz/eWr1aOgrC54MP3D/+0d15Z33yydYvfqGjUelAABpDsQM8rLo6OnlydOJEe8MG6SgIt88/N7fckn34YWv4cOu443RFhXQgANtGsQO8SG/cGH3+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d65mk5SQcS3L6dMne+ih2b59OeUKBB7FDgg+XV6+9RqLDRvcWbPM7NlmwQKdxx6W8D5XKeenP80OGOD0729at5aOA6BEKHZAiOg2bewTTlAnnGC++Wbz1KmR11+333nHouEFiLFt3bevPvTQmr32Mm3aSMcBUGoUOyCMdNu22SOPzB55pN68OTJ3rj13rv322xbfw/MtU1am+/Wz+ve3Dj5YtW6tlDJ57FwPwL8odkComYqKzODBmcGDVTZrv/9+ZN68yNy51ldfSedC8+y4ox4wQB98sPXTn6pYTDoNAHkUOwBKKaUiEad3b6d37/Tvfqc/+yyyYIG9YIG9cKGVzUonw/eYeFzvs4/u21cfcIDVubN0HADeQrED8ENm110zu+6aGTJEZbP2okX2ggWRd9+1Fi/W0sFCyyil99hDH3CA3n9/q2dPFucANIRiB6BhkYiz997O3nvXKaU3b7YWLYq8/779/vuUvBJwlbL23FP37q1797Z++lPFla0AmoFiB6BZTEWFs//+zv77K6WS8bj73nvmvffMBx+YRYu4QW2huImE06OH2727s9deTo8eyY4dpRMB8BmKHYCWSySsAw5QBxyglFKO437+ufngA/XBB2bRIvPpp+yB3Hyu1mbXXZ2f/MTZYw/nJz8xO+3EPSEA5INiByA/tm117qw6d1bHHKOUUum0WbbM/fhjtXix+egjs2QJdzP7LjceN126OF27Ot26uV27up07qwivwwAKhhcUAAUVj+sePewePbb+q+O4X36pli83y5aZZcvUsmXusmXh2TDPxOPuzjs7u+zi7rqru+uu7i67mA4dWJMDUDwUOwDFZNtWp06qUyd18MH1v5Fat06vW2etWqW//NJatar+f/rLL/3e9kxZme7USe28s955Z9Wxo+7USXfqpKuqUhs3SkcDECJ6/fr1OR+cTCZzPjaVx67ozGUucwM4d9Mmd/Vq9dVXZs0atXq1Wb1arVunUimzbp3auNELF+EapdzKStW2rdu2rdluO1NV5XboYNq3dzt0UFVVpqJim0d57v9n5jKXuYGey4odAG9o3dpq3Vrtvvs2/pPrmm++MevWbV65Um/erKqr9caNuqZGb96sNm/WtbUqnbZqa1Umo2prVTqt02nlOEop5Tgqm1Wuq+vvh2tZSmsViSjbVlor21ZlZbqsTJWVqVhMlZersjJVVqYrK1VFhaqsVK1b64oK1br1RttWbduaykrOogLwOIodAM+zLN2unW7Xzqmqyvkx8vkEzH1XAfgFnz4BAAACgmIHAAAQEBQ7AACAgKDYAQAABATFDgAAICAodgAAAAFBsQMAAAgIih0AAEBAUOwAAAACgmIHAAAQEBQ7AACAgKDYAQAABATFDgAAICAodgAAAAFBsQMAAAgIih0AAEBAUOwAAAACQhtjcj44lUrlfGwymWQuc5nLXOYyl7nMZW4B57JiBwAAEBAUOwAAgICg2AEAAAQExQ4AACAgKHYAAAABQbEDAAAICIodAABAQFDsAAAAAoJiBwAAEBAUOwAAgICg2AEAAATE/wdiEDEAH6MEZgAAAABJRU5ErkJggg=="
            
        AddCoords(mol) #este addCoords é melhor que o metodo default do rdkit para calcular posiçoes 2d, especialmente para info estequiometrica, macrociclos e cadeias longas
        img = Draw.MolToImage(mol=mol)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue())
        return img_str.decode("UTF-8")

    def save(self, path=None):
        """Guarda o html gerado interativamente, na localização e nome selecionados, abrindo janela de dialogo"""
        if path == None:
            path = asksaveasfilename(initialfile="resultados" + time.strftime("%Y_%m_%d_%H_%M_%S", time.gmtime()) + ".html", initialdir= "./",filetypes=[('HTML File', '*.html'), ('All Files', '*.*')])
        if path:
            shutil.copyfile(os.path.join(self.tempdir, "index.html"), path)

    def __del__(self):
        """pequeno override do __del__ so para apagar a pasta temporaria apos apagar o objeto com que esta relacionada"""
        shutil.rmtree(self.tempdir, ignore_errors=True)  # descomentar para apagar a pasta temporaria a cada execucao e poupar espaço





if __name__ == '__main__': #EXEMPLO
    disp = Mostrador() #cria um objeto
    data = pd.read_csv("./scripts/test_dataset.csv", sep=";") #le um ficheiro csv para pd.DataFrame
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
