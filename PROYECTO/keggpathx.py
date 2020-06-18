import pandas as pd
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#parsear XML
import csv   
from lxml import etree
#Parsear HTML
import requests
import urllib.request
from bs4 import BeautifulSoup
from urllib.request import urlopen
from bs4.builder import HTMLParserTreeBuilder

#Parsear HTML
import os


def dow(org):
    
    os.makedirs("Pathways/"+str(org))
    
    ORGANISM = org

    pathways = requests.get('http://rest.kegg.jp/list/pathway/' + ORGANISM)

    for line in pathways.text.split('\n'):
        pathwayid = line.split('\t')[0].replace('path:', '')
        kgml = requests.get('http://rest.kegg.jp/get/' + pathwayid + '/kgml')
        path = "Pathways/"+str(org)+"/"+pathwayid + '.xml'        
        with  open(path, 'w') as f:      # Descarga y guarda en la direcion definida como path
            f.write(kgml.text)
            f.close



def net(org, maxi):
    path=os.getcwd()
    os.chdir('Pathways/'+str(org))       #Cambia la ejecuicon del programa a otra carpeta. 
    ## Cargar archivos
    #print('Cual es el primer path por encima de 1000 para ',org)
    #maxi=input()
    ORGANISM = org

    pathways = requests.get('http://rest.kegg.jp/list/pathway/' + ORGANISM)  #lista de las vias metabolicas y otra informacion porpocionada por KEGG

    for line in pathways.text.split('\n'):         # Se recorre la lista propocionada por  KEGG donde se estraen los nombres de las vias de KEGG
        pathwayid = line.split('\t')[0].replace('path:', '') 
        A = pathwayid + ".xml"                  # Le ponemos la extencion .xml para poder cargar los archivos
        if A == str(org)+str(maxi)+'.xml':           
                break
        locals()[str(pathwayid)] = etree.parse(A) # Guarda los  XML con el nombre de la via asignado por KEGG
        locals()[str(pathwayid)] = locals()[str(pathwayid)].getroot()  # Para transformar todos lo XML en directorios de python
    
    ## PARA GUARDAR LOS NODOS
    NODESS=[]
    for line in pathways.text.split('\n'):        
        pathwayid = line.split('\t')[0].replace('path:', '') 
        
        if pathwayid == str(org)+str(maxi):           
                break
            
        ii=locals()[str(pathwayid)]
        
        locals()['NODES_'+str(pathwayid)]={}  #Este comando crea un directorio vacio.
        for i in range(0,len(ii)):
            if ii[i].tag=='entry':         
                A=ii[i].attrib            
                if A['type'] == 'gene' or 'group':
                    locals()['NODES_'+str(pathwayid)][A['id']]= A['name'] 
                 
                    NODESS.append( [ A['id']  ,  A['name']  ]  )


    ## PARA GUARDAR LOS VINCULOS
    EDGES=[]
    
    for line in pathways.text.split('\n'):        
        pathwayid = line.split('\t')[0].replace('path:', '') 
    
        if pathwayid == str(org)+str(maxi):          
                break
            
        jj=locals()[str(pathwayid)]
        
        
        #locals()['EDGES_'+str(pathwayid)]=[]
        #EDGES=locals()['EDGES_'+str(pathwayid)]
        NODES=locals()['NODES_'+str(pathwayid)]
        #print(NODES)
        for i in range(0,len(jj)):
            if jj[i].tag=='relation':         
                A=jj[i].attrib            
                if A['type'] == 'ECrel':    #relación enzima-enzima, que indica dos enzimas que catalizan pasos de reacción sucesivos
                    #print(A)
                    EDGES.append( [  NODES[ A['entry1'] ] , NODES[ A['entry2'] ]  ]  ) 
                    
    ## GRAFICAR GRAFO
    G = nx.DiGraph() #--->DiGrafo para hacer un grado dirigido    
    G.add_edges_from(EDGES)
    
    #volver carpeta
     
    os.chdir(path)       
    

    return G, NODESS , EDGES   


def net_label(org, maxi):
    path=os.getcwd()
    os.chdir('Pathways/'+str(org))       #Cambia la ejecuicon del programa a otra carpeta. 
    
    ORGANISM = org

    pathways = requests.get('http://rest.kegg.jp/list/pathway/' + ORGANISM)  #lista de las vias metabolicas y otra informacion porpocionada por KEGG

    for line in pathways.text.split('\n'):         # Se recorre la lista propocionada por  KEGG donde se estraen los nombres de las vias de KEGG
        pathwayid = line.split('\t')[0].replace('path:', '') 
        A = pathwayid + ".xml"                  # Le ponemos la extencion .xml para poder cargar los archivos
        if A == str(org)+str(maxi)+'.xml':           
                break
        locals()[str(pathwayid)] = etree.parse(A) # Guarda los  XML con el nombre de la via asignado por KEGG
        locals()[str(pathwayid)] = locals()[str(pathwayid)].getroot()  # Para transformar todos lo XML en directorios de python
    
    
    data=pd.read_csv('data_'+org+'.csv',index_col=0, sep=';')
    data=data[['n' ,'enzima', 'nob_Keeg' ,'path']]
    data=data.values
    
    for line in pathways.text.split('\n'):
    
        pathwayid = line.split('\t')[0].replace('path:', '') 
    
        if pathwayid == str(org)+str(maxi):           
                break
    
        locals()['NODES_'+str(pathwayid)]={}
    
        for i in data:
            if i[3] == pathwayid:                
                locals()['NODES_'+str(pathwayid)][i[0]]= i[1] 
                
    EDGES=[]

    for line in pathways.text.split('\n'):        
        pathwayid = line.split('\t')[0].replace('path:', '') 
    
        if pathwayid == str(org)+str(maxi):           
                break
            
        ii=locals()[str(pathwayid)]
    
        NODES=locals()['NODES_'+str(pathwayid)]
   
        for i in range(0,len(ii)):
            if ii[i].tag=='relation':         
                A=ii[i].attrib            
                if A['type'] == 'ECrel':    #relación enzima-enzima, que indica dos enzimas que catalizan pasos de reacción sucesivos
               
                    EDGES.append( [  NODES[ A['entry1'] ] , NODES[ A['entry2'] ]  ]  ) 

         
    ## GRAFICAR GRAFO
    G = nx.DiGraph() #--->DiGrafo para hacer un grado dirigido    
    G.add_edges_from(EDGES)
    
    #volver carpeta
     
    os.chdir(path)       
    

    return G, NODES , EDGES 


def group_betweenness(G):
    
    intr = nx.betweenness_centrality(G) 
    intr = pd.DataFrame([[key, intr[key]] for key in intr.keys()], columns=['Name', 'intermediacion'])
    intr = intr.sort_values(by='intermediacion', ascending=False)
    
    suma=0
    maxx=intr.values[0][1]
    for i in intr.values:
        suma=(maxx-i[1])+suma
        
    g=nx.number_of_nodes(G)
    
    return suma/(g-1)   


def group_cerc(G):
    
    cerc = nx.closeness_centrality(G)
    cerc = pd.DataFrame([[key, cerc[key]] for key in cerc.keys()], columns=['Name', 'Cercania'])                 
    cerc = cerc.sort_values(by='Cercania', ascending=False)
    
    suma=0
    maxx=cerc.values[0][1]
    for i in cerc.values:
        suma=(maxx-i[1])+suma
    
    g=nx.number_of_nodes(G)
    return suma/((g-1)*(g-2)/(2*g-3))
    
    
    
def enz_position_dfinout(G,name_enzima):

    d=G.degree(G) #Vector con el par ordenado de numero de aristas de entrada a cada uno de los nodos
    d_in=G.in_degree(G)
    d_out=G.out_degree(G)
    
    degre=[]
    for i in G.nodes:
        degre.append([i,d[i],str(d_in[i]), str(d_out[i]), d_in[i]-d_out[i] ])
    df=pd.DataFrame(degre,columns=['Nombre', 'grado' , 'grado_in', 'grado_out','dif_grado'])
    ac=df.sort_values(by='dif_grado', ascending=False)
   

    n=[]
    suma=0
    for i in range(0,len(ac)):
        suma=1+suma
        n.append(suma)

    ac['n']=n 

    #'enolase [EC:4.2.1.11]'
    num=ac[ac['Nombre'] == name_enzima ].values[0][4]

    
    return  num 
    
def enz_position_inter(G,name_enzima):
    
    intr=nx.betweenness_centrality(G) #dicionario en python
    intr = pd.DataFrame([[key, intr[key]] for key in intr.keys()], columns=['Nombre', 'intermediacion'])
    intr=intr.sort_values(by='intermediacion', ascending=False)

    n=[]
    suma=0
    for i in range(0,len(intr)):
        suma=1+suma
        n.append(suma)

    intr['n']=n 
    num=intr[intr['Nombre'] == name_enzima ].values[0][1]
    
    return num


def enzimas_out(EDGES,enzima):
    df=pd.DataFrame(EDGES, columns=['ENZIMA_ORIGEN', 'ENZIMA_DESTINO'])
    
    df=df[df['ENZIMA_ORIGEN'] == enzima ]
    df=df.drop_duplicates(subset = 'ENZIMA_DESTINO')
   
    df_out=df
    
    return df_out
    

def enzimas_in(EDGES,enzima):
    df=pd.DataFrame(EDGES, columns=['ENZIMA_ORIGEN', 'ENZIMA_DESTINO'])

    df=df[df['ENZIMA_DESTINO'] == enzima ]
    df=df.drop_duplicates(subset = 'ENZIMA_ORIGEN')
    df_in=df
    
    return df_in
    