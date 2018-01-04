# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 20:57:37 2017

@author: Erman
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

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
    cluster = nx.clustering(G)
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
    tris = nx.triangles(G)
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

#Inputs
N=23133
k=9
p=0.4

G_WS = nx.watts_strogatz_graph(N, k, p)


#Degree
PlotDegreeDistribution(G_WS, loglogplot=False)
print 'Average Network Degree', np.average(G_WS.degree().values())


#Hubs
percentage=0.05
print getHubsAndDegrees(G_WS,percentage)
print getHubs(G_WS,percentage)

#Clusters
PlotClusteringDistribution(G_WS, loglogplot=False)
print nx.average_clustering(G_WS)
print nx.clustering(G_WS)

#Triangles
G_WS_triangles=nx.triangles(G_WS)
print nx.triangles(G_WS)
PlotTriangleDistribution(G_WS, loglogplot=False)

#Edges and Nodes
print G_WS.number_of_edges()
print G_WS.number_of_nodes()

#Diameter
connComponents = nx.connected_component_subgraphs(G_WS)
for cc in connComponents:
    print cc.nodes()[0]
    print 'diameter of the connected component', nx.diameter(cc)
    print 'average shortest path length', nx.average_shortest_path_length(cc)


#Centralities
PlotDegreeCentralityDistribution(G_WS, loglogplot=False)
PlotClosenessCentralityDistribution(G_WS, loglogplot=False)
PlotBetweennessCentralityDistribution(G_WS, loglogplot=False)
PlotEigenvectorCentralityDistribution(G_WS, loglogplot=False)
PlotKatzCentralityDistribution(G_WS, loglogplot=False)

































