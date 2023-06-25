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
    l=sorted(l,key=lambda x:x[5])
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

ds=[]
#ds_cr=[]
with open("./sparse-0-0-ds-icp.csv") as f:
    for l in f:
        l=l.strip().split(",")
        l=list(map(convert,l))
        ds.append(l)



#with open("./sparse-0-cr-ds.csv") as f:
#    for l in f:
#        l=l.strip().split(",")
#        l=list(map(convert,l))
#        ds_cr.append(l)


a_size=[0]*2
anum=[0 for _ in range(2,40)]
ins=[0 for _ in range(25)]
ds=s_l(ds)
for l in ds:
    if l[5]==0.5:
        a_size[0]+=1
    if l[5]==4.5:
        a_size[1]+=1
    anum[l[0]-2]+=1
    if l[0]==7:
        ins[l[6]-1]+=1

print(a_size)
print([(i,anum[i-2]) for i in range(2,40)])
print([(i,ins[i-1]) for i in range(1,26)])
#ds_cr=s_l(ds_cr)

#for i in range(len(ds)):
#    if (ds[i][4]==1 and ds_cr[i][4]==0):
#        print("\n"+str(ds[i])+"\n"+str(ds_cr[i]))
'''
for i in range(len(v)):
    temp=[]
    for ll in All:
        temp.append(ll[i])
    
    check_SameCost(temp)
'''
