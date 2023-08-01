ori_data = {
    "sparse": {
        "van": [],
        "card": [],
        "ds": [],
        "card_ds": []
    },
    "dense": {
        "van": [],
        "card": [],
        "ds": [],
        "card_ds": []
    },
    "super-dense": {
        "van": [],
        "card": [],
        "ds": [],
        "card_ds": []
    }
}
our_data = {
    "sparse": {
        "van": [],
        "card": [],
        "ds": [],
        "card_ds": []
    },
    "dense": {
        "van": [],
        "card": [],
        "ds": [],
        "card_ds": []
    },
    "super-dense": {
        "van": [],
        "card": [],
        "ds": [],
        "card_ds": []
    }
}
algos=["van","card","ds","card_ds"]
ms=["sparse","dense","super-dense"]
f_algos={"van":"0-0","card":"card-0","ds":"0-ds","card_ds":"card-ds"}

def sortIns(l):
    l=sorted(l,key=lambda i:i[1])
    l=sorted(l,key=lambda i:i[0])
    return l


def conv(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

if __name__=="__main__":
    #read ori
    for algo in algos:
        with open(f"./CCBS_ORI/opt_{algo}.txt") as f:
            for l in f:
                for m in ms[::-1]:
                    if m in l:
                        l=l.replace(m,"").strip()
                        l=list(map(conv,l.split(",")[:-1]))
                        ori_data[m][algo].append(l)
                        break
    #read our
    for algo in algos:
        for m in ms:
            alg_name=f_algos[algo]
            with open(f"./CCBS_MY_ORI2/{m}-{alg_name}.csv") as f:
                for l in f:
                    l=l.strip()
                    l=list(map(conv,l.split(",")[:-1]))
                    our_data[m][algo].append(l)
    
    #sort based on agent number and ins
    for m,mv in ori_data.items():
        for algo,d in mv.items():
            ori_data[m][algo]=sortIns(ori_data[m][algo])
            our_data[m][algo]=sortIns(our_data[m][algo])

    #prelim check
    '''for m,mv in ori_data.items():
        for algo,d in mv.items():
            for i in range(len(ori_data[m][algo])):
                ori_l=ori_data[m][algo][i]
                our_l=our_data[m][algo][i]
                print(ori_l)
                print(our_l)
                assert ori_l[0]==our_l[0] and ori_l[1]==our_l[1]
                if ori_l[3]==1: 
                    assert our_l[3]==1 and our_l[5] == ori_l[5]
    '''

    #calc node expansion time
    ori_expand_total=our_expand_total=0
    for m,mv in ori_data.items():
        for algo,d in mv.items():
            for i in range(len(ori_data[m][algo])):
                ori_l=ori_data[m][algo][i]
                our_l=our_data[m][algo][i]

                if ori_l[2] >=29.99 and our_l[2]>=29.99:
                    ori_expand_total+=ori_l[-1]
                    our_expand_total+=our_l[-1]

    print("extra time for our node to expand: ",(ori_expand_total-our_expand_total)/ori_expand_total*100)


    our_lower_bound=[]
    ori_lower_bound=[]
    for m,mv in ori_data.items():
        for algo,d in mv.items():
            for i in range(len(ori_data[m][algo])):
                ori_l=ori_data[m][algo][i]
                our_l=our_data[m][algo][i]
                ori_lower_bound.append(ori_l[5])
                our_lower_bound.append(our_l[5])
    
    #[print(ori_lower_bound[i],our_lower_bound[i]) for i in range(len(our_lower_bound))]
    with open("lower Bound.csv",'w') as f:
        print(*ori_lower_bound,sep=',',file=f)
        print(*our_lower_bound,sep=',',file=f)

    


