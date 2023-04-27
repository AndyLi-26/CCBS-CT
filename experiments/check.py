def convert(x):
    try:
        x = int(x)
    except ValueError:
        try:
            x = float(x)
        except ValueError:
            pass
    return x
LL=[]
with open("./sparse-0-0-0.csv",'r') as f:
    for i,l in enumerate(f):
        l=l.strip()
        l=l.split(",")
        l=l[:-1]
        l=list(map(convert,l))
        if l[1]==4.5:
            LL.append(l)
    
    LL=sorted(LL,key=lambda x:x[2])
    [print(i) for i in LL]



