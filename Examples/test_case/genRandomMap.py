import random
import sys,math
#def Crossing(line1,line2,nodes): # return true if crossing, false if not
#    p1=nodes[line1[0]];p2=nodes[line1[1]]
#    p3=nodes[line2[0]];p4=nodes[line2[1]]
#    if p1[0]

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
