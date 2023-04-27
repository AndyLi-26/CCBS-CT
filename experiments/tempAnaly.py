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
with open("./batch1/super-dense.csv","r") as f:
    for i,l in enumerate(f):
        l=l.strip()
        l=l.split(",")
        l=l[:-1]
        l=list(map(convert,l))
        info.append(l)

info=sorted(info,key=lambda x:x[2])
#[print(i) for i in info]

n=len(info)/4
for r in [0.5,1.5,2.5,3.5,4.5]:
    print("r=",r)
    Vanilla=DS=CR=DS_CR=0
    V_succ=[0]*30
    D_succ=[0]*30
    C_succ=[0]*30
    DC_succ=[0]*30
    V_en=[0]*30
    D_en=[0]*30
    C_en=[0]*30
    DC_en=[0]*30
    V_c=[0]*30
    D_c=[0]*30
    C_c=[0]*30
    DC_c=[0]*30
    for i in info:
        if i[1]!=r: continue
        if not i[3] and not i[5]:
            Vanilla+=i[7]
            V_succ[i[0]]+=i[7]
            V_en[i[0]]+=i[11]
            V_c[i[0]]+=i[9]
        elif i[3] and not i[5]:
            DS+=i[7]
            D_succ[i[0]]+=i[7]
            D_en[i[0]]+=i[11]
            D_c[i[0]]+=i[9]
        elif not i[3] and i[5]:
            CR+=i[7]
            C_succ[i[0]]+=i[7]
            C_en[i[0]]+=i[11]
            C_c[i[0]]+=i[9]
        else:
            DS_CR+=i[7]
            DC_succ[i[0]]+=i[7]
            DC_en[i[0]]+=i[11]
            DC_c[i[0]]+=i[9]
    print("Vanilla:",Vanilla,Vanilla/n)
    [print(i,end=', ') for i in enumerate(V_succ)]
    [print(i,end=', ') for i in enumerate(V_en)]
    [print(i,end=', ') for i in enumerate(V_c)]
    print()
    print("DS:",DS,DS/n)
    [print(i,end=', ') for i in enumerate(D_succ)]
    [print(i,end=', ') for i in enumerate(D_en)]
    [print(i,end=', ') for i in enumerate(D_c)]
    print()
    print("CR:",CR,CR/n)
    [print(i,end=', ') for i in enumerate(C_succ)]
    [print(i,end=', ') for i in enumerate(C_en)]
    [print(i,end=', ') for i in enumerate(C_c)]
    print()
    print("DS+CR:",DS_CR,DS_CR/n)
    [print(i,end=', ') for i in enumerate(DC_succ)]
    [print(i,end=', ') for i in enumerate(DC_en)]
    [print(i,end=', ') for i in enumerate(DC_c)]
    print()
    with open("r_"+str(int(r)),'w') as f:
        [print(i,end=', ',file=f) for i in (V_succ)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (V_en)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (V_c)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (D_succ)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (D_en)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (D_c)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (C_succ)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (C_en)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (C_c)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (DC_succ)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (DC_en)]
        print("",file=f)
        [print(i,end=', ',file=f) for i in (DC_c)]
        print("",file=f)




inslist=[0]*26
alist=[0]*30
sizelist=[0]*5
temp1=[0]*26
temp2=[0]*5
temp3=0
for l in info:
    inslist[l[2]]+=1
    alist[l[0]]+=1
    sizelist[int(l[1])]+=1
    if l[0]==5 and l[1]==4.5:
        temp1[l[2]]+=1
    if l[0]==5:
        temp2[int(l[1])]+=1

    if l[8]==-1:
        temp3+=1

'''
[print(i,end=', ') for i in enumerate(inslist)]
print()
[print(i,end=', ') for i in enumerate(alist)]
print()
print(sizelist)

[print(i,end=' ') for i in enumerate(temp1)]
print()
[print(i,end=' ') for i in enumerate(temp2)]
print()
print(temp3)
'''
#print("--------------")
#temp=sorted(info,key=lambda x:x[2])
#newinfo=[]
#for l in temp:
    #print(l[2]-1,Alllist[l[2]-1])
#    if Alllist[l[2]-1]:
#        newinfo.append(l)
#print("--------------")
#[print(i) for i in newinfo]


#print(alst)
