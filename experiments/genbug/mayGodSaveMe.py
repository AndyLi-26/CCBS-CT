#from itertools import combinations
import subprocess,os
import random
import math
import shutil
nodes=[]
dist=lambda x,y:math.sqrt((nodes[x][0]-nodes[y][0])**2+(nodes[x][1]-nodes[y][1])**2)
N=10
P=5
processPoll=[]
lst=list(range(170))
idx_taken=[]
class idx:
    def __init__(self):
        self.idx=0
        self.idx_take=[]

    def __call__(self,using):
        while 1:
            self.idx=(self.idx+1)%10000000
            if (not self.idx in self.idx_take) and (not self.idx in using):
                return self.idx
    def add_taken(self,i):
        self.idx_take.append(i)

idx_iter=idx()

def readMap():
    with open("../../Instances/roadmaps/sparse/map.graph") as f:
        node_num=int(f.readline())
        for i in range(node_num):
            l=f.readline()
            nodes.append(list(map(float,l.split(" "))))

def overlapping(l):
    for i in range(N):
        for j in range(i+1,N):
            if dist(l[i],l[j])<=1.5:
                return True
    return False

def ok(l1,l2):
    for i in range(N):
        if l1[i]==l2[i]:
            return False
    if overlapping(l1): return False
    if overlapping(l2): return False
    return True



def getNewTasks():
    while 1:
        l1=random.sample(lst,k=N)
        l2=random.sample(lst,k=N)
        if ok(l1,l2):
            return l1,l2

def getNewTaskFile(idx):
    l1,l2=getNewTasks() 
    if os.path.exists(f"./{idx}"):
        shutil.rmtree(f"./{idx}")
    os.makedirs(f"./{idx}", exist_ok=True)
    with open(f"./{idx}/task.task",'w') as f:
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
    cmd1=["../../release/CCBS","-m","../../Instances/roadmaps/sparse/map.graph",
         "-t",f"./{idx}/task.task",
        "--HI_h","0",
         "-o", f"{idx}/result",
         "-a","0.5",
         "--agent_num",str(N),
         "--timelimit","5",
         "--extra_info",str(idx),
         "--DS","--CT"]

    cmd2=["../../release/CCBS","-m","../../Instances/roadmaps/sparse/map.graph",
         "-t",f"./{idx}/task.task",
        "--HI_h","0",
         "-o", f"{idx}/result",
         "-a","0.5",
         "--agent_num",str(N),
         "--timelimit","5",
         "--extra_info",str(idx),
         "--DS"]
    return [cmd1,cmd2]

def checkResult(idx):
    with open(f"./{idx}/result") as f:
        res=[]
        for l in f:
            l=l.strip().split(",")
            res.append(l)
        
        assert len(res)==2,"err file len at "+str(idx)
        if res[0][8]=='1' and res[1][8]=='1':
            if float(res[0][10])-float(res[1][10])>1e-5:
                idx_iter.add_taken(idx)
                with open("./bug_file","a") as f:
                    f.write(f'{idx}: DSCT: {res[1][-3]}, DS: {res[0][-3]}\n')
                return

    shutil.rmtree(f"./{idx}")


        
if __name__=="__main__":
    readMap()
    processPool=[-1]*P
    pending=[[] for i in range(P)]
    idxRec=[-1]*P
    while 1:
        for i in range(P):
            if processPool[i]==-1 or processPool[i].poll() is not None:
                if len(pending[i])==0:
                    if idxRec[i]!=-1:
                        print(idxRec[i])
                        print(processPool[i])
                        checkResult(idxRec[i])
                    newIDX=idx_iter(idxRec)
                    getNewTaskFile(newIDX)
                    pending[i]=getNewCmds(newIDX)
                    idxRec[i]=newIDX
                cmd=pending[i].pop()
                print(subprocess.list2cmdline(cmd))
                processPool[i]=subprocess.Popen(cmd)
                print("just printing",processPool[i])
    #temp_combination=list(combinations(lst,4))


