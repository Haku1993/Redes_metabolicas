#Parsear HTML
import os
import requests

def way(org):
    
    os.makedirs("Pathways/"+str(org))
    
    ORGANISM = org

    pathways = requests.get('http://rest.kegg.jp/list/pathway/' + ORGANISM)

    for line in pathways.text.split('\n'):
        pathwayid = line.split('\t')[0].replace('path:', '')
        kgml = requests.get('http://rest.kegg.jp/get/' + pathwayid + '/kgml')
        path = "Pathways/"+str(org)+"/"+pathwayid + '.xml'        
        with  open(path, 'w') as f:
            f.write(kgml.text)
            f.close
