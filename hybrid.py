# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 09:38:48 2017

@author: Erman
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random as rd
import operator

#Node Degrees

def PlotDegreeDistribution(G, loglogplot=False):
    degrees = G.degree()
    values = sorted(degrees.values())
    hist = [degrees.values().count(x) for x in values]
    #plot degree distribution
    plt.figure()
    plt.grid(True)
    if not loglogplot:
        plt.plot(values, hist, 'bv-')
    else:
        plt.loglog(values, hist, 'bv-')
    plt.legend(['Degree - k'])
    plt.xlabel('Degree')
    plt.ylabel('Number Of Nodes')
    plt.title('Degree Distribution of the Network')
    #plt.savefig(os.path.join(pathout, figName))
    #plt.close()
    plt.show()

def getHubsAndDegrees(G,difPercentage):
    revDeg = []
    for k,v in G.degree().items():
        revDeg.append((v,k))
    degSorted = sorted(revDeg, reverse=True)
    
    hubsDegrees=[]
    valPrevious=degSorted[0][0]
    for tup in degSorted:
        value=tup[0]
        key=tup[1]
#        print value, valPrevious,(float(valPrevious-value)/float(valPrevious))
        if (float(valPrevious-value)/float(valPrevious))<difPercentage:
            hubsDegrees.append((key,value))
        else:
            break
        valPrevious=value
    return hubsDegrees

def getHubs(G,difPercentage):
    revDeg = []
    for k,v in G.degree().items():
        revDeg.append((v,k))
    degSorted = sorted(revDeg, reverse=True)
    
    hubs=[]
    valPrevious=degSorted[0][0]
    for tup in degSorted:
        value=tup[0]
        key=tup[1]
#        print value, valPrevious,(float(valPrevious-value)/float(valPrevious))
        if (float(valPrevious-value)/float(valPrevious))<difPercentage:
            hubs.append(key)
        else:
            break
        valPrevious=value
    return hubs
#Clusters
def PlotClusteringDistribution(G, loglogplot=False):
    cluster = nx.clustering(G_original)
    values = sorted(cluster.values())
    hist = [cluster.values().count(x) for x in values]
    #plot degree distribution
    plt.figure()
    plt.grid(True)
    if not loglogplot:
        plt.plot(values, hist, 'bv-')
    else:
        plt.loglog(values, hist, 'bv-')
    plt.legend(['Clustering - k'])
    plt.xlabel('clustering')
    plt.ylabel('Number Of Nodes')
    plt.title('Clustering Distribution of the Network')
    #plt.savefig(os.path.join(pathout, figName))
    #plt.close()
    plt.show()


    
#Triangles
def PlotTriangleDistribution(G, loglogplot=False):
    tris = nx.triangles(G_original)
    values = sorted(tris.values())
    hist = [tris.values().count(x) for x in values]
    #plot degree distribution
    plt.figure()
    plt.grid(True)
    if not loglogplot:
        plt.plot(values, hist, 'bv-')
    else:
        plt.loglog(values, hist, 'bv-')
    plt.legend(['Triangles - k'])
    plt.xlabel('Triangles')
    plt.ylabel('Number Of Nodes')
    plt.title('Triangle Distribution of the Network')
    #plt.savefig(os.path.join(pathout, figName))
    #plt.close()
    plt.show()
    
#Degree Centrality
def PlotDegreeCentralityDistribution(G, loglogplot=False):
    degCent = nx.degree_centrality(G)
    values = sorted(degCent.values())
    hist = [degCent.values().count(x) for x in values]
    #plot degree distribution
    plt.figure()
    plt.grid(True)
    if not loglogplot:
        plt.plot(values, hist, 'bv-')
    else:
        plt.loglog(values, hist, 'bv-')
    plt.legend(['Degree Centrality'])
    plt.xlabel('Degree Centralities')
    plt.ylabel('Number Of Nodes')
    plt.title('Degree Centrality of the Network')
    #plt.savefig(os.path.join(pathout, figName))
    #plt.close()
    plt.show()
    
#Closeness Centrality
    
def PlotClosenessCentralityDistribution(G, loglogplot=False):
    clsCent = nx.closeness_centrality(G)
    values = sorted(clsCent.values())
    hist = [clsCent.values().count(x) for x in values]
    #plot degree distribution
    plt.figure()
    plt.grid(True)
    if not loglogplot:
        plt.plot(values, hist, 'bv-')
    else:
        plt.loglog(values, hist, 'bv-')
    plt.legend(['Closeness Centrality'])
    plt.xlabel('Closeness Centralities')
    plt.ylabel('Number Of Nodes')
    plt.title('Closeness Centrality of the Network')
    #plt.savefig(os.path.join(pathout, figName))
    #plt.close()
    plt.show()    
    
#Closeness Centrality
    
def PlotBetweennessCentralityDistribution(G, loglogplot=False):
    btwCent = nx.betweenness_centrality(G)
    values = sorted(btwCent.values())
    hist = [btwCent.values().count(x) for x in values]
    #plot degree distribution
    plt.figure()
    plt.grid(True)
    if not loglogplot:
        plt.plot(values, hist, 'bv-')
    else:
        plt.loglog(values, hist, 'bv-')
    plt.legend(['Betweenness Centrality'])
    plt.xlabel('Betweenness Centralities')
    plt.ylabel('Number Of Nodes')
    plt.title('Betweenness Centrality of the Network')
    #plt.savefig(os.path.join(pathout, figName))
    #plt.close()
    plt.show()  
    
#EigenVector Centrality
    
def PlotEigenvectorCentralityDistribution(G, loglogplot=False):
    eigCent = nx.eigenvector_centrality(G)
    values = sorted(eigCent.values())
    hist = [eigCent.values().count(x) for x in values]
    #plot degree distribution
    plt.figure()
    plt.grid(True)
    if not loglogplot:
        plt.plot(values, hist, 'bv-')
    else:
        plt.loglog(values, hist, 'bv-')
    plt.legend(['Eigenvector Centrality'])
    plt.xlabel('Eigenvector Centralities')
    plt.ylabel('Number Of Nodes')
    plt.title('Eigenvector Centrality of the Network')
    #plt.savefig(os.path.join(pathout, figName))
    #plt.close()
    plt.show()  
    
#Katz Centrality
    
def PlotKatzCentralityDistribution(G, loglogplot=False):
    katzCent = nx.katz_centrality(G)
    values = sorted(katzCent .values())
    hist = [katzCent .values().count(x) for x in values]
    #plot degree distribution
    plt.figure()
    plt.grid(True)
    if not loglogplot:
        plt.plot(values, hist, 'bv-')
    else:
        plt.loglog(values, hist, 'bv-')
    plt.legend(['Katz Centrality'])
    plt.xlabel('Katz Centralities')
    plt.ylabel('Number Of Nodes')
    plt.title('Katz Centrality of the Network')
    #plt.savefig(os.path.join(pathout, figName))
    #plt.close()
    plt.show()     


#Read Nodes in

def readNodes(nodeList,fileName):
    with open(fileName, "r") as ins:
        for line in ins:
            line=line.split()
            node1=int(line[0])
            node2=int(line[1])
            
            nodeList.append(node1)
            
            nodeList.append(node2)
       
        nodeList=list(set(nodeList))
    return nodeList   

#Add Node          
def addRandom(G,node,k):
    nodes=G.nodes()
    N=len(nodes)
    p=float(k/(float(N+1)))
    
    G.add_node(node)
       
    for item in nodes:
        number=rd.random()
        if number<p:
            G.add_edge(node,item)
            
def addPreferential(G,node,k):
    nodes=G.nodes() 
    dSum=G.number_of_edges()       
    
    G.add_node(node)
    
    for item in nodes:
        number=rd.random()
        p=k*(float(G.degree(item))/(float(dSum)+1.0))
        if number<p:
            G.add_edge(node,item)

def generateHybrid(G,nodeList,kAverage,pSelection):            
    ii=0            
    for item in nodeList:
        number=rd.random()
        if number < pSelection:
           addRandom(G,item,kAverage) 
#           print "random"
        else:
           addPreferential(G,item,kAverage)
#           print "prefer"
        ii+=1
        print ii
            
#Main 
N=23133
#Average expected k
kAverage=8.0
#Node Lists
nodeList=[]
#File Names
fileName='CA-CondMat.txt'
#Graph
G_hybrid=nx.Graph()
#Probabilities
pSelection=0.5


#Get Node Lists            
nodeList=readNodes(nodeList,fileName)
#Create the graph
generateHybrid(G_hybrid,nodeList,kAverage,pSelection)
print nx.info(G_hybrid)

#Hubs
percentage=0.03
print getHubsAndDegrees(G_hybrid,percentage)
print getHubs(G_hybrid,percentage)


#Degree
PlotDegreeDistribution(G_hybrid, loglogplot=False)
print 'Average Network Degree', np.average(G_hybrid.degree().values())


#Clusters
PlotClusteringDistribution(G_hybrid, loglogplot=False)
print nx.average_clustering(G_hybrid)
print nx.clustering(G_hybrid)

#Triangles
G_hybrid_triangles=nx.triangles(G_hybrid)
print nx.triangles(G_hybrid)
PlotTriangleDistribution(G_hybrid, loglogplot=False)

#Edges and Nodes
print G_hybrid.number_of_edges()
print G_hybrid.number_of_nodes()


#Centralities
PlotDegreeCentralityDistribution(G_hybrid, loglogplot=False)
PlotClosenessCentralityDistribution(G_hybrid, loglogplot=False)
PlotBetweennessCentralityDistribution(G_hybrid, loglogplot=False)
PlotEigenvectorCentralityDistribution(G_hybrid, loglogplot=False)
PlotKatzCentralityDistribution(G_hybrid, loglogplot=False)


######################################################################

#Average Clustering
print nx.average_clustering(G_hybrid)

#Triangle Number
tris = nx.triangles(G_hybrid)
print sum(tris.values())/3

#Diameters
diameters=[]
connComponents = nx.connected_component_subgraphs(G_hybrid)
for cc in connComponents:
     diameters.append( nx.diameter(cc))
     
diameters.sort(reverse=True)
print diameters[0]


#Connected Components
print nx.number_connected_components(G_hybrid)


#Largest Component Analysis
Gcc=sorted(nx.connected_component_subgraphs(G_hybrid), key = len, reverse=True)
G0=Gcc[0]

print "Number of Nodes: ",G0.number_of_nodes()
print "Number of Edges: ",G0.number_of_edges()

print "Average path length",nx.average_shortest_path_length(G0)


#First Two Connected
Gcc=sorted(nx.connected_component_subgraphs(G_hybrid), key = len, reverse=True)
G0=Gcc[0]
G1=Gcc[1]

print "Ratio of nodes in the second largest to the largest: ", float(Gcc[1].number_of_nodes())/Gcc[0].number_of_nodes()
print "Ratio of edgees in the second largest to the largest: ", float(Gcc[1].number_of_edges())/Gcc[0].number_of_edges()


print "Ratio of nodes in the  largest to the whole: ", float(Gcc[0].number_of_nodes())/G_hybrid.number_of_nodes()
print "Ratio of edgees in the  largest to the whole: ", float(Gcc[0].number_of_edges())/G_hybrid.number_of_edges()


#First Five  Centralities
number=5
degCent = nx.degree_centrality(G_hybrid)
sorted_Cent = sorted(degCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_Cent[ii]
    print tupl



clsCent = nx.closeness_centrality(G_hybrid)
sorted_Cls = sorted(clsCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_Cls[ii]
    print tupl


btwCent = nx.betweenness_centrality(G_hybrid)
sorted_btw = sorted(btwCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_btw[ii]
    print tupl


eigCent = nx.eigenvector_centrality(G_hybrid)
sorted_eig = sorted(eigCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_eig[ii]
    print tupl


katzCent = nx.katz_centrality_numpy(G_hybrid)
sorted_katz = sorted(katzCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_katz[ii]
    print tupl


nx.draw(G_hybrid, nodelist=clsCent.keys(), node_size=[10000*v for v in clsCent.values()])
plt.show()

nx.draw(G_hybrid, nodelist=btwCent.keys(), node_size=[100000*v for v in btwCent.values()])
plt.show()


nx.draw(G_hybrid, nodelist=eigCent.keys(), node_size=[100000*v for v in eigCent.values()])
plt.show()

nx.draw(G_hybrid, nodelist=katzCent.keys(), node_size=[50000*v for v in katzCent.values()])
plt.show()
    
    
    
#Page Rank
pr = nx.pagerank(G_hybrid)
number=5

sorted_Pr = sorted(pr.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_Pr[ii]
    print tupl

#Hits and Authorities
h,a=nx.hits(G_hybrid)

number=5
print "hubs nodes"

sorted_h = sorted(h.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_h[ii]
    print tupl

number=5
print "authority nodes"

sorted_a = sorted(a.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_a[ii]
    print tupl


    
#Cliques
ratioDenominator=1000
smallestSize=nx.number_of_nodes(G_hybrid)/ratioDenominator

communities=list(nx.k_clique_communities(G_hybrid, smallestSize))

for community in communities:
        print list(community)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

