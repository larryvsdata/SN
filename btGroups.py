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
    plt.legend(['Degree '])
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
    plt.legend(['Clustering '])
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
    plt.legend(['Triangles '])
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
#        ii+=1
#        print ii
            
def generateHybridPercentage(G,nodeList,kAverage,pSelection,startPercent,endPercent):            
    start=int(startPercent*len(nodeList))      
    end=int(endPercent*len(nodeList))      
    for index in range(start,end):
        item=nodeList[index]
        number=rd.random()
        if number < pSelection:
           addRandom(G,item,kAverage) 
#           print "random"
        else:
           addPreferential(G,item,kAverage)
           
def generateHInGBg(G,nodeList,kAverage,pSelection,groupPercent,intervalP,intervalS,percentageP,percentageS):            
    groupStartIndex=int(len(nodeList)*groupPercent)
    
    for index in range(len(nodeList)):
        item=nodeList[index]
        number=rd.random()
        if number < pSelection:
           addRandom(G,item,kAverage) 
#           print "random"
        else:
           addPreferential(G,item,kAverage)           
           
        if index>groupStartIndex and index%intervalP==0:
            
            hubs=getHubs(G,percentageP)
            print hubs
            for hub in hubs:
                print "Hub ",hub
                group=G.neighbors(hub)
                connectList2(G,group)
                
        if index>groupStartIndex and index%intervalS==0:
            print "Inter hubs connection"
            connectHubs(G,percentageS)
                
        print "Index ",index  
                
                
#Connect Nodes in the list
#
#def connectList(G,nodeList,pConnection):
#    
#    for ii in range(len(nodeList)-1):
#        for jj in range(ii+1,len(nodeList)):
#            number=rd.random()
#            if number < pConnection:
#                G.add_edge(nodeList[ii],nodeList[jj])
                
def connectList2(G,nodeList):
    index1=rd.randint(0,len(nodeList)-1)
    index2=rd.randint(0,len(nodeList)-1)
    
    if index1!=index2:
        G.add_edge(nodeList[index1],nodeList[index2])

                
#def connectTwoLists(G,list1,list2,pConnection):
#    
#    for item1 in list1:
#        for item2 in list2:
#            number=rd.random()
#            if number < pConnection:
#                G.add_edge(item1,item2)
#                
                
def connectTwoLists2(G,list1,list2):
    number1=rd.randint(0,len(list1)-1)
    number2=rd.randint(0,len(list2)-1)

    G.add_edge(list1[number1],list2[number2])     


def  connectHubs(G,percentage):
     hubList=getHubs(G,percentage) 
     if len(hubList)==2:
         list1=G.neighbors(hubList[0])
         list2=G.neighbors(hubList[1])
         connectTwoLists2(G,list1,list2)
     elif len(hubList)>2:
         maximum=20
         if len(hubList)>2:
             hubList=hubList[0:maximum]
         for ii in range(len(hubList)-1):
             for jj in range(ii+1,len(hubList)):
                 print ii,jj,len(hubList)
                 list1=G.neighbors(hubList[ii])
                 list2=G.neighbors(hubList[jj])
                 connectTwoLists2(G,list1,list2)
            


#Main 
N=23133
#Average expected k
kAverage=8.0
#Node Lists
nodeList=[]
#File Names
fileName='CA-CondMat.txt'
#Graph
G_hybrid_btG=nx.Graph()
#Probabilities
pSelection=0.5
# When do add-ons start
groupPercent=0.5
#Intervals
intervalPrimary=10
intervalSecondary=10
#Percentage to decide about Hubs
percentageP=0.07
percentageS=0.15   
#Get Node Lists        
nodeList=readNodes(nodeList,fileName)
#Create the graph
generateHInGBg(G_hybrid_btG,nodeList,kAverage,pSelection,groupPercent,intervalPrimary,intervalSecondary,percentageP,percentageS)

print nx.info(G_hybrid_btG)

#Hubs
percentage=0.07
print getHubsAndDegrees(G_hybrid_btG,percentage)
print getHubs(G_hybrid_btG,percentage)


#Degree
PlotDegreeDistribution(G_hybrid_btG, loglogplot=False)
print 'Average Network Degree', np.average(G_hybrid_btG.degree().values())


#Clusters
PlotClusteringDistribution(G_hybrid_btG, loglogplot=False)
print nx.average_clustering(G_hybrid_btG)
print nx.clustering(G_hybrid_btG)

#Triangles
G_hybrid_triangles=nx.triangles(G_hybrid_btG)
print nx.triangles(G_hybrid_btG)
PlotTriangleDistribution(G_hybrid_btG, loglogplot=False)

#Edges and Nodes
print G_hybrid_btG.number_of_edges()
print G_hybrid_btG.number_of_nodes()


#Centralities
PlotDegreeCentralityDistribution(G_hybrid_btG, loglogplot=False)
PlotClosenessCentralityDistribution(G_hybrid_btG, loglogplot=False)
PlotBetweennessCentralityDistribution(G_hybrid_btG, loglogplot=False)
PlotEigenvectorCentralityDistribution(G_hybrid_btG, loglogplot=False)
PlotKatzCentralityDistribution(G_hybrid_btG, loglogplot=False)


######################################################################

#Average Clustering
print nx.average_clustering(G_hybrid_btG)

#Triangle Number
tris = nx.triangles(G_hybrid_btG)
print sum(tris.values())/3

#Diameters
diameters=[]
connComponents = nx.connected_component_subgraphs(G_hybrid_btG)
for cc in connComponents:
     diameters.append( nx.diameter(cc))
     
diameters.sort(reverse=True)
print diameters[0]


#Connected Components
print nx.number_connected_components(G_hybrid_btG)


#Largest Component Analysis
Gcc=sorted(nx.connected_component_subgraphs(G_hybrid_btG), key = len, reverse=True)
G0=Gcc[0]

print "Number of Nodes: ",G0.number_of_nodes()
print "Number of Edges: ",G0.number_of_edges()

print "Average path length",nx.average_shortest_path_length(G0)


#First Two Connected
Gcc=sorted(nx.connected_component_subgraphs(G_hybrid_btG), key = len, reverse=True)
G0=Gcc[0]
G1=Gcc[1]

print "Ratio of nodes in the second largest to the largest: ", float(Gcc[1].number_of_nodes())/Gcc[0].number_of_nodes()
print "Ratio of edgees in the second largest to the largest: ", float(Gcc[1].number_of_edges())/Gcc[0].number_of_edges()


print "Ratio of nodes in the  largest to the whole: ", float(Gcc[0].number_of_nodes())/G_hybrid_btG.number_of_nodes()
print "Ratio of edgees in the  largest to the whole: ", float(Gcc[0].number_of_edges())/G_hybrid_btG.number_of_edges()


#First Five  Centralities
number=5
degCent = nx.degree_centrality(G_hybrid_btG)
sorted_Cent = sorted(degCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_Cent[ii]
    print tupl



clsCent = nx.closeness_centrality(G_hybrid_btG)
sorted_Cls = sorted(clsCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_Cls[ii]
    print tupl


btwCent = nx.betweenness_centrality(G_hybrid_btG)
sorted_btw = sorted(btwCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_btw[ii]
    print tupl


eigCent = nx.eigenvector_centrality(G_hybrid_btG)
sorted_eig = sorted(eigCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_eig[ii]
    print tupl


katzCent = nx.katz_centrality_numpy(G_hybrid_btG)
sorted_katz = sorted(katzCent.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_katz[ii]
    print tupl


nx.draw(G_hybrid_btG, nodelist=clsCent.keys(), node_size=[10000*v for v in clsCent.values()])
plt.show()

nx.draw(G_hybrid_btG, nodelist=btwCent.keys(), node_size=[100000*v for v in btwCent.values()])
plt.show()


nx.draw(G_hybrid_btG, nodelist=eigCent.keys(), node_size=[100000*v for v in eigCent.values()])
plt.show()

nx.draw(G_hybrid_btG, nodelist=katzCent.keys(), node_size=[50000*v for v in katzCent.values()])
plt.show()
    
    
    
#Page Rank
pr = nx.pagerank(G_hybrid_btG)
number=5

sorted_Pr = sorted(pr.items(), key=operator.itemgetter(1), reverse=True)
for ii in range(number):
    tupl=sorted_Pr[ii]
    print tupl

#Hits and Authorities
h,a=nx.hits(G_hybrid_btG)

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
smallestSize=nx.number_of_nodes(G_hybrid_btG)/ratioDenominator

communities=list(nx.k_clique_communities(G_hybrid_btG, smallestSize))

for community in communities:
        print list(community)
    







