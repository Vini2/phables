# Source: https://github.com/algbio/MFD-ILP

import itertools
import math as m
import networkx as nx
# import matplotlib.pyplot as plt


def read_input(graphfile):
    
    trip_data = open(graphfile,'r').read().split('\n')
    i = 0
    listOfGraphs = {}
    k = 0
    
    while(True):
        if "#" in trip_data[i]:
            i = i+1
            N = int(trip_data[i])
            edges = list()
            while(True):
                i = i+1;
                if "#" in trip_data[i]:
                    break;
                if "" == trip_data[i]:
                    break;
                line = trip_data[i].split(" ")
                edges.append((line[0],line[1],line[2],line[3]))
                if i >= len(trip_data)-1:
                    break;
            G = {'Nodes':N,'list of edges': edges}
            listOfGraphs[k] = G
            k +=1 
        if i >= len(trip_data)-1:
                    break; 

    return listOfGraphs;


def flowMultipleDecomposition(data,K):
    
   
    # calculate the minimal flow decomposition based on such graph
    f_low = data['flows_low']
    f_up = data['flows_up']
    V = data['vertices']
    E = data['edges']
    W = data['maxFlow']
    S = data['sources']
    D = data['targets']
    AD_in = data['adj_in']
    AD_out = data['adj_out']
    
    # import library
    from docplex.mp.model import Model

    # model name
    Flow = Model(name = 'Flow Decomposition')

    # define variables
    x = {(i,j,k): Flow.binary_var(name="x_{0}_{1}_{2}".format(i,j,k)) for (i,j) in E for k in range(0,K)}
    w = {(k): Flow.integer_var(name="w_{0}_".format(k)) for k in range(0,K)}
    z = {(i,j,k): Flow.integer_var(name="z_{0}_{1}_{2}".format(i,j,k)) for (i,j) in E for k in range(0,K)}

    # flow conservation
    for k in range(0,K):
        for i in V:
            if i in S:
                Flow.add_constraint(sum(x[i,j,k] for j in AD_out[i])  == 1) 
            if i in D:
                 Flow.add_constraint(sum(x[j,i,k] for j in AD_in[i]) == 1)
            if i not in S and i not in D:
                Flow.add_constraint(sum(x[i,j,k] for j in AD_out[i]) - sum(x[j,i,k] for j in AD_in[i]) == 0)

    
    # flow balance
    for (i,j) in E:
        Flow.add_constraint(f_up[i,j] >= sum(z[i,j,k] for k in range(0,K)))
        Flow.add_constraint(f_low[i,j] <= sum(z[i,j,k] for k in range(0,K)))
        
    
    # linearization
    for (i,j) in E:
        for k in range(0,K):
            Flow.add_constraint(z[i,j,k] <= W*x[i,j,k])
            Flow.add_constraint(w[k] - (1 - x[i,j,k])*W <= z[i,j,k])
            Flow.add_constraint(z[i,j,k] <= w[k])
    
    # objective function
    Flow.minimize(1)

    Flow.parameters.threads = 1
    Flow.parameters.mip.display = 0
    Flow.parameters.timelimit = 3600
    
    sol = Flow.solve(log_output=False)
    print('status',Flow.solve_details.status)

    if Flow.solve_details.status == 'integer infeasible' or Flow.solve_details.status == 'time limit exceeded, no integer solution':
        data['messsage'] = "infeasible"
        return data

    elif Flow.solve_details.status == 'time limit exceeded':
        print('time limit')
  
    print('FD time',Flow.solve_details.time)

    # calculate averages
    x_sol = {(i,j,k):Flow.solution.get_value(x[(i,j,k)]) for (i,j) in E for k in range(0,K)}
    w_sol = {(k):Flow.solution.get_value(w[(k)]) for (k) in range(0,K)}
    paths =[list() for i in range(0,K)]
    for k in range(0,K):
        for (i,j) in E:
            if x_sol[i,j,k] == 1:
                paths[k].append((i,j))
               
    data['weights'] = w_sol
    data['solution'] = paths
    data['message'] = 'solved'
    data['runtime'] = Flow.solve_details.time
    
    return data;
    
    

def FD_Algorithm(data):
    
    listOfEdges = data['edges']
    solutionMap = data['graph']
    solutionSet = 0
    Kmin = data['minK']
    solutionWeights = 0

    for i in range(1,len(listOfEdges)+1):
        data = flowMultipleDecomposition(data,i)
        if data['message'] == "solved":
            solutionSet = data['solution']
            solutionWeights = data['weights']
            break;

    # Get solution paths and weights
    solution_paths = {}
    
    if solutionSet != 0:
        for i in range(0,len(solutionSet)):
            solution_paths[i] = {"weight": solutionWeights[i], "path": solutionSet[i]}
            # print("W:",solutionWeights[i], solutionSet[i])
    
    return data, solution_paths


def SolveInstances(Graphs,outfile,recfile):
    
    fp = open(outfile,'w+')
    fc = open(recfile,'w+')
    
    for i in range(0,1): 
        
        f_low = {}
        f_up = {}
        Edges = set()
        V = set()
        listOfEdges = Graphs[i]['list of edges']
        for k in range(0,len(listOfEdges)):
            (a,b,c,d) = (listOfEdges[k])
            Edges.add((a,b))
            V.add(a)
            V.add(b)
            f_low[a,b] = int(float(c))
            f_up[a,b] = int(float(d))
            
        # creation of graphs
        G = nx.DiGraph()
        G.add_edges_from(Edges,weights = f_low)
        G.add_nodes_from(V)
        
        # creation of adjacent matrix
        AD_in = {}
        AD_out = {}
        
        for v in V:
            setAdj = set()
            for (i,j) in list(G.out_edges(v)):
                if i != v:
                    setAdj.add(i)
                if j != v:
                    setAdj.add(j)
            
            AD_out[v] = list(setAdj)
            
            setAdj = set()
            for (i,j) in list(G.in_edges(v)):
                if i != v:
                    setAdj.add(i)
                if j != v:
                    setAdj.add(j)
            
            AD_in[v] = list(setAdj)
            
        
        # calculating source, sinks and max flows
        S = [x for x in G.nodes() if G.out_degree(x)>=1 and G.in_degree(x)==0]
        D = [x for x in G.nodes() if G.out_degree(x)==0 and G.in_degree(x)>=1]
        maxW = max(f_up.values())
        
        
        # definition of data
        
        data = {'edges' : Edges,
                'flows_low' : f_low,
                'flows_up': f_up,
                'vertices' : V,
                'graph' : G,
                'Kmax' : len(Edges),
                'weights' : {},
                'sources': S,
                'targets': D,
                'message': {},
                'solution': 0,
                'maxFlow': maxW,
                'adj_in': AD_in,
                'adj_out': AD_out,
                'minK': 2,
                'runtime': 0,
        }
        
        data, solution_paths = FD_Algorithm(data)
    
    return solution_paths


