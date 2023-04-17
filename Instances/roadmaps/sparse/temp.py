with open("map.graph",'r') as f:
    n=int(f.readline())
    nodes=[]
    for _ in range(n):
        nodes.append(f.readline())

    n=int(f.readline())
    edges=[]
    for l in range(n):
        l=f.readline().split(' ')
        l=list(map(int,l))
        if l[0]==85: l[0]=120
        if l[1]==85: l[1]=120
        if l in edges:
            print(l)
        else:
            edges.append(l)


with open("new_map.graph","w") as f:
    print(len(nodes),file=f)
    [print(i,end='',file=f) for i in nodes]
    print(len(edges),file=f)
    [print(*i,file=f) for i in edges]
