import sys
if __name__=="__main__":
    with open("map.graph") as f:
        node_num=int(f.readline().strip())
        nodes=[]
        for i in range(node_num):
            info=f.readline().strip()
            info=tuple(map(float,info.split(' ')))
            nodes.append(info)
        edge_num=int(f.readline().strip())
        edges=[]
        for i in range(edge_num):
            info=f.readline().strip()
            info=tuple(map(int,info.split(' ')))
            edges.append(info)

    [print(i) for i in nodes]
    [print(i) for i in edges]



    with open("nodes.csv","w") as f:
        [print(*i,sep=",",file=f) for i in nodes]

    
    with open("edges.csv","w") as f:
        [print(*i,sep=",",file=f) for i in edges]
