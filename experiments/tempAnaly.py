def convert(x):
    try:
        x = int(x)
    except ValueError:
        try:
            x = float(x)
        except ValueError:
            pass
    return x
info=[]
with open("sparse.csv","r") as f:
    for i,l in enumerate(f):
        l=l.strip()
        l=l.split(",")
        l=l[:-1]
        l=list(map(convert,l))
        info.append(l)

info=sorted(info,key=lambda x:x[2])
[print(i) for i in info]

Vanilla=DS=CR=DS_CR=0
for i in info:
    if not i[3] and not i[5]:
        Vanilla+=i[7]
    elif i[3] and not i[5]:
        DS+=i[7]
    elif not i[3] and i[5]:
        CR+=i[7]
    else:
        DS_CR+=i[7]

print("Vanilla:",Vanilla,Vanilla/25)
print("DS:",DS,DS/25)
print("CR:",CR,CR/25)
print("DS+CR:",DS_CR,DS_CR/25)

inslist=[0]*25
alist=[0]*6
sizelist=[0]*5
for l in info:
    inslist[l[2]-1]+=1
    alist[l[0]-5]+=1
    sizelist[int(l[1])]+=1

[print(i,end=', ') for i in enumerate(inslist)]
print()
print(alist)
print(sizelist)
#print("--------------")
temp=sorted(info,key=lambda x:x[2])
newinfo=[]
for l in temp:
    #print(l[2]-1,Alllist[l[2]-1])
    if Alllist[l[2]-1]:
        newinfo.append(l)
#print("--------------")
#[print(i) for i in newinfo]


print(alst)
