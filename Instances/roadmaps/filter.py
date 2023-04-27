import random
import math
N=30

def getNewTasks(Pool):
    temp=[random.sample(Pool,k=N),random.sample(Pool,k=N)]
    return [(temp[0][i],temp[1][i]) for i in range(N)]


if __name__=="__main__":
    dist=lambda x,y:math.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2)
    
    for name in ["dense"]:
        with open(name+"/map.graph") as f:
            nodes=[]
            node_num=int(f.readline())
            for i in range(node_num):
                l=f.readline()
                nodes.append(list(map(float,l.split(" "))))
            #[print(i) for i in nodes]
            #print("----------------------------------------------------")
        
        Pool=[i for i in range(len(nodes)) if i not in [355,356,522,549,683,759,760]]

        for ind in range(1,26):
            taskName=name+"/"+str(ind)+"task.task"
            '''with open(taskName) as f:
                pairs=[]
                for i in f:
                    pairs.append(list(map(int,i.split(" "))))
                #print(pairs)
            '''

            while 1:
                pairs=getNewTasks(Pool)
                for i in range(N):
                    #print(i)
                    for j in range(N):
                        if i==j:
                            continue
                        dis1=dist(nodes[pairs[i][0]],nodes[pairs[j][0]]) 
                        dis2=dist(nodes[pairs[i][1]],nodes[pairs[j][1]])
                        #print(dis1,dis2)
                        if dis1 <10:
                            break
                        if dis2 <10:
                            break
                    else:
                        continue
                    break
                else:
                    break
            print(ind," pass")
            print(sorted(pairs,key=lambda x:x[0]))
            with open(taskName,'w') as f:
                print(N,file=f)
                [print(*i,file=f) for i in pairs]
        
