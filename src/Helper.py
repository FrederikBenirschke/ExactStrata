#Helper functions
import sage

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method

from sage.modules.free_module import FreeModule  # pylint: disable=import-error
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from exactstrata.profile import *
from exactstrata.exactboundarystratum import *
from exactstrata.iteratedblowup import *
from exactstrata.divtautpolynomial import *
from sage.calculus.var import var






def marked_to_component(X):
    '''Dictionary listing the marked points of a generalised stratum. The marked points in the ambient are enumerated by [1...n].
    The value for a kay is a tuple (i,j) pointing to the j-th marked point on the i-th component.''' 
    return { i+j : (i,j) for i,sig in enumerate(X.sig_list()) for j,_ in enumerate(sig.sig)}


            

    
    
def path_between_markings(ELG, start, end):
    '''Returns a path between two marked points on the dual graph (implemented via BFS on the dual graph).'''
    G=ELG.LG
    return bfs(G,start,end)

 
def bfs(LG, start_leg,goal_leg):
    '''# Breadth-first search to find a path between two vertices/marked points in the dual graph.
    Assumes that both vertices are in the same connected component of the graph.'''
    visited = [] # List with entries (vertex,leg) to keep track of visited nodes.
    path_dict= {}
    queue = []     # Initialize a queue
    visited.append([LG.vertex(start_leg),start_leg, 0 ])

    queue.append(LG.vertex(start_leg))

    while queue:
        s = queue.pop( 0 ) 

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
                        if V[ 0 ]==position:
                            path+=[V[ 1 ],V[ 2 ]]
                            break
                    position = LG.vertex(path[- 1 ])

                path.append(start_leg)
                path.reverse()
                return path
