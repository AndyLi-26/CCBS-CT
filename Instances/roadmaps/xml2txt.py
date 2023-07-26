import sys,os,shutil
def remove_rubbish(s):
    x=s.index("<agent")
    s=s[x:]
    s=s.strip()
    return s
    
def parserT(info):
    b=info.split('"')
    t=(int(b[1]),int(b[3]))
    return t
    
def parserP(i,info):
    x=info.index("<path")
    info=info[x:]
    info=list(map(lambda x:x.strip(), info.split(">")))[:-2]
    dur=info[0].split('"')[1]
    section=info[1:]
    with open("paths/{}.csv".format(i),'w') as f:
        i=[0,0,0,0,dur]
        print(*i,sep=',',file=f)
        for i in section:
            i=i.split('"')[3::2]
            print(*i,sep=',',file=f)
    return dur
    #with open("paths/{}.csv".format(i),'w'):
    
if __name__=="__main__":
    for i in range(1,26):
        with open("super-dense/ori_set_xml/"+str(i)+"_task.xml",'r') as f:
            s=f.read().strip()
        s=remove_rubbish(s)
        all=list(map(lambda x:x.strip(),s.split("<log>")))
        #tasks
        tasks=list(map(lambda x:x.strip(),all[0].split("<agent")))[1:]
        with open("super-dense/ori_set/"+str(i)+"_task.task",'w') as f:
            print(len(tasks),file=f)
            for info in tasks:
                print(*parserT(info),sep=' ',file=f)
                
