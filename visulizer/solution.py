import sys,os,shutil
paths=[]
with open(sys.argv[1]) as f:
    for l in f:
        l=l.strip()
        info=l.split(",")
        if len(info)==2:
            i=int(info[0])
            paths.append([])
            paths[i].append("0,0,0,0,"+info[1]+",")
        else:
            paths[i].append(l)

[print(i) for i in paths]

isExist = os.path.exists('paths')
print(isExist)
if isExist:
    shutil.rmtree("paths")
os.mkdir("./paths")

for i,p in enumerate(paths):
    print("---------",i,"---------")
    with open("paths/{}.csv".format(i),'w') as f:
        [print(l,file=f) for l in p]
