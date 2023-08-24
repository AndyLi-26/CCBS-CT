import json
def check_SameCost(infos):
    temp=[i[10] for i in infos if i[8]]
    if len(temp)<2: return True
    a=temp[0]
    for i in temp[1:]:
        if abs(a-i)>0.00001:
            print("---")
            [print(i) for i in infos]
            [print(i) for i in temp]
            print("---")

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

with open("./config.json","r") as f:
    config=json.loads(f.read())
    for k,v in config.items():
        config[k]=v.split(" ")
    print(config)
output_address="{}-{}-{}-{}-{}.csv"
algo_name="{}-{}-{}-{}"
for m in config['m']:
    all_data[m]=dict()
    for es in config['es']:
        for ct in config['ct']:
            for ds in config['ds']:
                for icp in config['icp']:
                    es_tag="es" if es=='1' else '0'
                    if ct=='0': ct_tag="0" 
                    elif ct=='1': ct_tag="ct"
                    elif ct=="2": ct_tag="ct_abs"
                    ds_tag="ds" if ds=='1' else '0'
                    icp_tag="icp" if icp=='1' else '0'
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
for a in config['a']:
    for r in config['r']:
        for i in config['i']:
            idx_list.append((int(a),float(r),int(i)))



for m,l in all_data.items():
    for idx in idx_list:
        temp=[]
        for k,v in l.items():
            if idx in v:
                temp.append(v[idx])
        #print("----")
        #[print(i) for i in temp]
        #if idx[1]==4.5: continue
        #if m=="sparse" and idx==(5,4.5,70): continue
        check_SameCost(temp)


[[print(m, k,len(v)) for k,v in i.items()] for m,i in all_data.items()]





