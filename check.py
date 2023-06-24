def check_SameCost(infos):
    temp=[i[10] for i in infos if i[8]]
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
    ins_idx=l[6]
    if ins_idx not in ins:
        ins[ins_idx]=[]
    ins[ins_idx].append(l)

for k in ins:
    print("("+str(k)+":",str(len(ins[k]))+")")
    [print((i[10] if i[8] else ""),i[1:5]+i[8:]) for i in ins[k]]
    #for i in ins[k]:
    #    print(i[4],": ")
    #    print(i[1:4]+i[7:])






