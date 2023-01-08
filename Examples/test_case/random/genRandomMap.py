import random
import sys,math
#def Crossing(line1,line2,nodes): # return true if crossing, false if not
#    p1=nodes[line1[0]];p2=nodes[line1[1]]
#    p3=nodes[line2[0]];p4=nodes[line2[1]]
#    if p1[0]

s_node='''    <node id="n{}">
      <data key="key0">{},{}</data>
    </node>'''

def exportNode(nodes,f):
    for i,(x,y) in enumerate(nodes):
        print(s_node.format(i,x,y))
        print(s_node.format(i,x,y),file=f)

s_edge='''    <edge id="e{}" source="n{}" target="n{}">
      <data key="key1">1</data>	
    </edge>'''

def exportEdge(edges,f):
    for i,(e1,e2) in enumerate(edges):
        print(s_edge.format(i,e1,e2))
        print(s_edge.format(i,e2,e1))
        print(s_edge.format(i,e1,e2),file=f)
        print(s_edge.format(i,e2,e1),file=f)

def exportMap(nodes,edges):
    s_beg='''<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="key0" for="node" attr.name="coords" attr.type="string" />
  <key id="key1" for="edge" attr.name="weight" attr.type="double" />
  <graph id="G" edgedefault="directed" parse.nodeids="free" parse.edgeids="canonical" parse.order="nodesfirst">'''

    s_end='''  </graph>
</graphml>
    '''
    with open('map.xml','w') as f:
        print(s_beg,file=f)
        print(s_beg)
        exportNode(nodes,f)
        exportEdge(edges,f)
        print(s_end,file=f)
        print(s_end)

def exportTask(tasks):
    s_beg='''<?xml version="1.0" ?>
<root>'''
    s_tasks='''   <agent start_id="{}" goal_id="{}"/>'''
    s_end='''</root>'''
    with open('tasks.xml','w') as f:
        print(s_beg,file=f)
        print(s_beg)
        for s,e in tasks:
            print(s_tasks.format(s,e),file=f)
            print(s_tasks.format(s,e))
        print(s_end,file=f)
        print(s_end)

if __name__=="__main__":
    nodes=[(0,0)]
    r=5
    level=1
    min_dis=4*r
    max_num=6
    N=int(sys.argv[1])
    holes_left=max_num
    cur_level=[]
    dis=lambda n1,n2 : math.sqrt((n1[0]-n2[0])**2+(n1[1]-n2[1])**2)
    #gen nodes:
    while len(nodes)<N:
        circumf=min_dis*level*8
        for _ in range(100):
            pos=random.random()*circumf
            idx=pos//(min_dis*level*2)
            match int(idx):
                case 0:
                    x=-level*min_dis
                    y=-level*min_dis+pos%(min_dis*level*2)
                case 1:
                    x=-level*min_dis+pos%(min_dis*level*2)
                    y=level*min_dis
                case 2:
                    x=level*min_dis
                    y=level*min_dis-pos%(min_dis*level*2)
                case 3:
                    x=level*min_dis-pos%(min_dis*level*2)
                    y=-level*min_dis

            for n in cur_level:
               if dis(n,(x,y))<min_dis:
                   break
            else:
                cur_level.append((x,y))
                nodes.append((x,y))
                break
        else:
            level+=1
            cur_level=[]
    minx=miny=10
    for x,y in nodes:
        minx=min(x,minx)
        miny=min(y,miny)
    final=[(x+abs(minx)*1.1+r,y+abs(miny)*1.1+r) for (x,y) in nodes]
    random.shuffle(final)
    print(len(final))
    [print(*n,sep=',') for n in final]
    with open("nodes.csv",'w') as f:
        [print(*n,sep=',',file=f) for n in final]

    #gen edge
    edge=[]
    covered=[0]
    uncovered=list(range(1,N))
    for _ in range(1,N):
        min_dis=10000000000
        min_idx=(-1,)
        for i in covered:
            for j in uncovered:
                temp=dis(final[i],final[j])
                if temp<min_dis:
                    min_dis=temp
                    min_idx=(i,j)
        uncovered.remove(min_idx[1])
        covered.append(min_idx[1])
        edge.append(min_idx)
    
    #add 50% of random
    while len(edge)<int(1.5*N):
        temp=tuple(random.sample(range(N),2))
        if temp in edge or (temp[1],temp[0]) in edge:
            continue
        else:
            edge.append(temp)

    print(*edge)
    with open('edges.csv','w') as f:
        [print(*n,sep=',',file=f) for n in edge]

    #gen tasks
    while 1:
        A=random.sample(range(N),N)
        B=random.sample(range(N),N)
        tasks=[]
        for i in range(N):
            if A[i]==B[i]:
                break
            else:
                tasks.append((A[i],B[i]))
        else:
            break
    print(*tasks)
    exportMap(final,edge)
    exportTask(tasks)
