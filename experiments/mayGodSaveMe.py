#from itertools import combinations
import random
import math
import os
nodes=[]
dist=lambda x,y:math.sqrt((nodes[x][0]-nodes[y][0])**2+(nodes[x][1]-nodes[y][1])**2)
N=5
processPoll=[]
lst=list(range(170))
idx_taken=[]
def readMap():
    with open("../Instances/roadmaps/sparse/map.graph") as f:
        node_num=int(f.readline())
        for i in range(node_num):
            l=f.readline()
            nodes.append(list(map(float,l.split(" "))))
def overlapping(l):
    for i in range(N):
        for j in range(i+1,N):
            if dist(i,j)<=1:
                return True
    return False

def ok(l1,l2):
    for i in range(N):
        if l1[i]==l2[i]:
            return False
    if overlapping(l1): return False
    if overlapping(l2): return False



def getNewTasks():
    while 1:
        l1=random.sample(lst,k=N)
        l2=random.sample(lst,k=N)
        if ok(l1,l2):
            return l1,l2

def regen(n):
   l1,l2=getNewTasks() 
    with open("./for_bugs/"+n+"_task.task",'w') as f:
        print(N,file=f)
        [print(l1[i],l2[i],file=f) for i in range(N)]

def addcmd(cms):
    if (len(processPool)>=9):
        finish = False
        while not finish:
            time.sleep(1)
            for p in range(0,len(processPool)):
                if p >= len(processPool):
                    break
                if processPool[p].poll() is not None:
                    processPool.pop(p)
                    finish = True
                    p-=1
    else:
        for p in range(0,len(processPool)):
                if p >= len(processPool):
                    break
                if processPool[p].poll() is not None:
                    processPool.pop(p)
                    finish = True
                    p-=1
    
    try:
        processPool.append(subprocess.Popen(cmd))
    except:
        print(len(processPool))

def getNewCmds(idx):
    cmd1=[exe,"-m","../Instances/roadmaps/sparse/map.graph",
         "-t","./for_bugs/"+idx+"_task.task",
        "--HI_h","0",
         "-o", str(idx)+txt,
         "-a",0.5,
         "--agent_num",5,
         "--timelimit","120",
         "--extra_info",idx,
         "--DS","--CT"]

    cmd2=[exe,"-m","../Instances/roadmaps/sparse/map.graph",
         "-t","./for_bugs/"+idx+"_task.task",
        "--HI_h","0",
         "-o", "./for_bugs/"+str(idx)+".result",
         "-a",0.5,
         "--agent_num",5,
         "--timelimit","120",
         "--extra_info",idx,
         "--DS"]
    return [cmd1,cmd2]

def checkResult(idx):
    with open("./for_bugs/"+str(idx)+".result") as f:
        res=[]
        for l in f:
            l=l.strip().split(",")
            res.append(l)
        
        assert(len(l)==2,"err file len at "+idx)
        if res[0][8]=='1' and res[1][8]=='1':
            if res[0][10]!=res[1][10]:
                idx_taken.append(idx)
                with open("bug file:","a") as f:
                    f.write(f'{idx}: {res[0][-2]}, {res[1][-2]}\n')
        #remove files


        
if __name__=="__main__":
    readMap()
    idx=0
    processPool=[-1]*30
    pending=[[] for i in range(30)]
    idxRec=[-1]*30
    while 1:
        for i in range(30):
            if processPool[i]==-1 or processPool[i].poll() is not None:
                if len(pending[i])==0:
                    if idxRec[i]!=-1:
                        checkResult(idxRec[i])

                    pending[i]=getNewCmds(idx)
                    idx+=1
                cmd=pending.pop()
                processPool.append(subprocess.Popen(cmd))
    #temp_combination=list(combinations(lst,4))


