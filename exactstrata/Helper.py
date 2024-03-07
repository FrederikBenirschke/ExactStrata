#Helper functions


# Dictionary listing the marked points of a generalised stratum     
def marked_to_component(X):
    return { i+j : (i,j) for i,sig in enumerate(X.sig_list()) for j,_ in enumerate(sig.sig)}


            

    
    
def path_between_markings(ELG, start, end):
    G=ELG.LG
    return bfs(G,start,end)

# Breadth-first search to find a path between two vertices
# Only works if both vertices are in the same connected component of the graph  
def bfs(LG, start_leg,goal_leg):
    visited = [] # List of (vertex,leg) to keep track of visited nodes.
    path_dict= {}
    queue = []     # Initialize a queue
    visited.append([LG.vertex(start_leg),start_leg,_sage_const_0 ])

    queue.append(LG.vertex(start_leg))

    while queue:
        s = queue.pop(_sage_const_0 ) 

    #print (s, end = " ") 
        neighbors = []
        for i,j in  LG.edgesatvertex(s):
            if LG.vertex(i)== s:
                neighbors.append([LG.vertex(j),j,i])
            if LG.vertex(j)==s:
                neighbors.append([LG.vertex(i),i,j])

        for vertex,leg_i,leg_o in neighbors:
            if vertex not in [visit_vertex for visit_vertex,_,_ in visited]:
                visited.append([vertex,leg_i,leg_o])
                queue.append(vertex)
            if vertex==LG.vertex(goal_leg):

                # Now need to compute the path from goal_leg to start_leg
                position=vertex

                path=[goal_leg]
                while(position!=LG.vertex(start_leg)):


                    for V in visited:
                        if V[_sage_const_0 ]==position:
                            path+=[V[_sage_const_1 ],V[_sage_const_2 ]]
                            break
                    position = LG.vertex(path[-_sage_const_1 ])

                path.append(start_leg)
                path.reverse()
                return path
