def convert(x):
    try:
        x = int(x)
    except ValueError:
        x = float(x)
    return x
info=[]
with open("sparse.csv","r") as f:
    for i,l in enumerate(f):
        l=l.strip()
        l=l.split(",")
        l=l[:-1]
        l=list(map(convert,l))
        info.append(l)

Alllist=[4]*25


for l in info:
    Alllist[l[2]-1]-=1
[print(i,item) for i, item in enumerate(Alllist)]

print("--------------")
temp=sorted(info,key=lambda x:x[2])
newinfo=[]
for l in temp:
    print(l[2]-1,Alllist[l[2]-1])
    if Alllist[l[2]-1]:
        newinfo.append(l)
print("--------------")
[print(i) for i in newinfo]
