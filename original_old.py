# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 18:12:46 2017

@author: Erman
"""

import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
from collections import OrderedDict



def getReverseDict( Dict):
    reverseDict={}
    for word in Dict:
        wValue=Dict[word]
        if wValue not in reverseDict:
            reverseDict[wValue]=[word]
        else:
            reverseDict[wValue].append(word)
    return reverseDict
    
def incrementToDict(key,Dictionary):
    if key in Dictionary:
        Dictionary[key]+=1
    else:
        Dictionary[key]=1

    return Dictionary  
    
def getClusterClasses(clusters): 
    clusterValues=clusters.values()
    valueDict={}
    for value in clusterValues:
        valueDict=incrementToDict(value,valueDict)
    return valueDict
        
def getClusters(valueDict,maxNum):
     valueKeys=valueDict.keys()
     valueKeys.sort(reverse=True)
     valueKeys=valueKeys[0:maxNum]
     numbers=[]
     for val in valueKeys:
         numbers.append(valueDict[val])
     return [numbers,valueKeys]

#     plt.bar( numbers,valueKeys, align='center', alpha=0.5)
#
#     plt.ylabel('Usage')
#     plt.title('Programming language usage')
# 
#     plt.show()
     
def getNodeCountWMaxCluster(valueDict):  
    return valueDict[1.0]
    
         
     
    

#Read in    
G_original=nx.read_edgelist('CA-CondMat.txt',create_using=nx.Graph(),nodetype=int)


#Clusters
Avclustering_original=nx.average_clustering(G_original)
clustering_original=nx.clustering(G_original)

vDict=getClusterClasses(clustering_original)
maxNum=20
getClusters(vDict,maxNum)


#Triangles
G_original_triangles=nx.triangles(G_original)

#Info
print nx.info(G_original)




#degree_sequence=sorted(nx.degree(G_original).values(),reverse=True)