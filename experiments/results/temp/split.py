def check_SameCost(infos):
    temp=[i[9] for i in infos if i[7]]
    if len(temp)<2: return True
    a=temp[0]
    for i in temp[1:]:
        if abs(a-i)>0.00001:
            [print(i) for i in infos]
            print("---")
            [print(i) for i in temp]
            assert False
    



def s_l(l):
    l=sorted(l,key=lambda x:x[2])
    l=sorted(l,key=lambda x:x[0])
    l=sorted(l,key=lambda x:x[1])
    return l

def convert(x):
    try:
        x = int(x)
    except ValueError:
        try:
            x = float(x)
        except ValueError:
            pass
    return x

v=[]
ds=[]
cr=[]
ds_cr=[]
with open("./sparse-0-0-0.csv") as f:
    for l in f:
        l=l.strip().split(",")
        l=list(map(convert,l))
        v.append(l)
with open("./sparse-0-0-ds.csv") as f:
    for l in f:
        l=l.strip().split(",")
        l=list(map(convert,l))
        ds.append(l)
with open("./sparse-0-cr-0.csv") as f:
    for l in f:
        l=l.strip().split(",")
        l=list(map(convert,l))
        cr.append(l)
with open("./sparse-0-cr-ds.csv") as f:
    for l in f:
        l=l.strip().split(",")
        l=list(map(convert,l))
        ds_cr.append(l)


v,ds,cr,ds_cr=map(s_l,[v,ds,cr,ds_cr])
with open("./sparse-0-0-0.csv","w") as f:
    [print(*l,sep=',',file=f) for l in v]
with open("./sparse-0-0-ds.csv","w") as f:
    [print(*l,sep=',',file=f) for l in ds]
with open("./sparse-0-cr-0.csv","w") as f:
    [print(*l,sep=',',file=f) for l in cr]
with open("./sparse-0-cr-ds.csv","w") as f:
    [print(*l,sep=',',file=f) for l in ds_cr]

All=[v,ds,cr,ds_cr]
print(All)
assert len(v)==len(ds)==len(cr)==len(ds_cr),"len is not the same"
for i in range(len(v)):
    temp=[]
    for ll in All:
        temp.append(ll[i])
    
    check_SameCost(temp)






