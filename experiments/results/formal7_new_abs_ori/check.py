import json
def check_SameCost(infos):
    temp=[i[10] for i in infos if i[8]]
    if len(temp)<2: return 0
    a=temp[0]
    for i in temp[1:]:
        if abs(a-i)>0.00001:
            print("---")
            [print(i) for i in infos]
            [print(i) for i in temp]
            print("---")
            return 1
    return 0

def s_l(l):
    l=sorted(l,key=lambda x:x[6])
    l=sorted(l,key=lambda x:x[0])
    l=sorted(l,key=lambda x:x[5])
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

all_data=dict()

output_address="{}-{}-{}-{}-{}.csv"
algo_name="{}-{}-{}-{}"
algo=["","ds","ds ct1","ds ct2","ds icp","ds icp ct2"]
for m in ["sparse","dense","super-dense"]:
    all_data[m]=dict()
    for alg in algo:
        es_tag="es" if "es" in alg else '0'
        if "ct1" in alg: ct_tag="ct"
        elif "ct2" in alg: ct_tag="ct_abs"
        else: ct_tag="0"
        ds_tag="ds" if "ds" in alg else '0'
        icp_tag="icp" if "icp" in alg else '0'
        
        fName=output_address.format(m,es_tag,ct_tag,ds_tag,icp_tag)
        with open(fName) as f:
            temp=dict()
            for l in f:
                l=l.strip().split(",")
                l=list(map(convert,l))
                idx=(l[0],l[5],l[6])
                temp[idx]=l
            algoName=algo_name.format(es_tag,ct_tag,ds_tag,icp_tag)
            all_data[m][algoName]=temp
idx_list=[]
for a in range(3,41):
    for r in ["0.353553385"]:
        for i in range(1,26):
            idx_list.append((int(a),float(r),int(i)))

err=0
for m,l in all_data.items():
    for idx in idx_list:
        temp=[]
        for k,v in l.items():
            if idx in v:
                temp.append(v[idx])
            else:
                print("missing: ",m,k,idx)
        #print("----")
        #[print(i) for i in temp]
        #if idx[1]==4.5: continue
        #if m=="sparse" and idx==(5,4.5,70): continue
        err+=check_SameCost(temp)



[[print(m, k,len(v)) for k,v in i.items()] for m,i in all_data.items()]
print("total_err: ",err)




