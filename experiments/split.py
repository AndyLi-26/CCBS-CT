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

m="./batch1/sparse.csv"
v=[]
ds=[]
cr=[]
ds_cr=[]
with open(m,'r') as f:
    for i,l in enumerate(f):
        l=l.strip()
        l=l.split(",")
        l=l[:-1]
        l=list(map(convert,l))
        if not l[3] and not l[5]:
            v.append(l)
        elif l[3] and not l[5]:
            ds.append(l)
        elif not l[3] and l[5]:
            cr.append(l)
        else:
            ds_cr.append(l)

v,ds,cr,ds_cr=map(s_l,[v,ds,cr,ds_cr])
All=[v,ds,cr,ds_cr]
assert len(v)==len(ds)==len(cr)==len(ds_cr),"len is not the same"
for i in range(len(v)):
    temp=[]
    for ll in All:
        temp.append(ll[i])
    
    check_SameCost(temp)






