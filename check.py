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

result=[]
with open("./test_result.txt") as f:
    for l in f:
        l=l.strip().split(",")
        l=list(map(convert,l))
        result.append(l)

ins=dict()
for l in result:
    ins_idx=l[3]
    if ins_idx not in ins:
        ins[ins_idx]=[]
    ins[ins_idx].append(l)

for k in ins:
    print("("+str(k)+":",str(len(ins[k]))+")")
    [print((i[7] if i[5] else ""),i[1:2]+i[5:]) for i in ins[k]]
    #for i in ins[k]:
    #    print(i[4],": ")
    #    print(i[1:4]+i[7:])






