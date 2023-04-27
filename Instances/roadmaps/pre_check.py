import heapq
import math
def dis(nodes,n1,n2):
    return math.sqrt((nodes[n1][0]-nodes[n2][0])**2+(nodes[n1][1]-nodes[n2][1])**2)

def read_map(f_n):
    N_node=int(f_n.readline())
    Nodes=[]
    print(N_node)
    for _ in range(N_node):
        temp=f_n.readline()
        #print(temp)
        n=tuple(map(float,temp.split(" ")))
        Nodes.append(n)

    N_edge=int(f_n.readline())
    Edges=[[-1 for _ in range(N_node)] for _ in range(N_node)]
    for _ in range(N_edge):
        e=tuple(map(int,f_n.readline().split(" ")))
        Edges[e[0]][e[1]]=dis(Nodes,*e)
    
    #print(Edges)
    return Nodes,Edges



def dijkstra(graph, start):
    # Initialize distances and visited vertices
    print(len(graph))
    distances = [float('inf')]*len(graph)
    distances[start] = 0
    visited = set()

    # Initialize heap
    heap = [(0, start)]

    while heap:
        # Pop vertex with smallest distance from heap
        current_distance, current_vertex = heapq.heappop(heap)

        # Skip vertex if it has already been visited
        if current_vertex in visited:
            continue

        # Visit vertex
        visited.add(current_vertex)

        # Update distances of adjacent vertices
        for neighbor, weight in enumerate(graph[current_vertex]):
            if weight==-1: continue
            distance = current_distance + weight
            if distance < distances[neighbor]:
                distances[neighbor] = distance
                heapq.heappush(heap, (distance, neighbor))

    #print(distances)
    return distances


if __name__=="__main__":
    with open("./dense/map.graph","r") as f:
        N,G=read_map(f)
        P=dijkstra(G,0)
        for i,d in enumerate(P):
            if d==float('inf'):
                print(i,d,N[i])

