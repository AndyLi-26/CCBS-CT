import random
import math
N=90

dist=lambda x,y:math.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2)

def genList(Pool,nodes):
    retval=[]
    random.shuffle(Pool)
    while 1:
        if len(retval)>=N:
            return retval
        new_val=Pool.pop()
        for i in retval:
            if dist(nodes[i],nodes[new_val])<9.1:
                break
        else:
            retval.append(new_val)

def getNewTasks(Pool):
    temp=[random.sample(Pool,k=N),random.sample(Pool,k=N)]
    return [(temp[0][i],temp[1][i]) for i in range(N)]


if __name__=="__main__":
    noPool=[[85],[355,356,522,549,683,759,760],[]] 
    for m,name in enumerate(["sparse","dense","super-dense"]):
        with open(name+"/map.graph") as f:
            nodes=[]
            node_num=int(f.readline())
            for i in range(node_num):
                l=f.readline()
                nodes.append(list(map(float,l.split(" "))))
        
        Pool=[i for i in range(len(nodes)) if i not in noPool[m]]

        for ind in range(1,101):
            taskName=name+"/large_agent/"+str(ind)+"_task.task"
            starts=genList(Pool[:],nodes[:])
            goals=genList(Pool[:],nodes[:])
            while 1:
                for i in range(N):
                    if starts[i]==goals[i] or dist(nodes[starts[i]],nodes[goals[i]]) < 10:
                        random.shuffle(starts)
                        random.shuffle(goals)
                        break
                else:
                    print(ind," pass")
                    with open(taskName,'w') as f:
                        print(N,file=f)
                        [print(starts[i],goals[i],file=f) for i in range(N)]
                    break
            
