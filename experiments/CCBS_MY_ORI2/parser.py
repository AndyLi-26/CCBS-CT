data = {
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
def conv(s):
    try:
        return int(s)
    except ValueError:
        return float(s)
f_algos={"van":"0-0","card":"card-0","ds":"0-ds","card_ds":"card-ds"}

if __name__=="__main__":
    for algo in algos:
        for m in ms:
            alg_name=f_algos[algo]
            with open(f"{m}-{alg_name}.csv") as f:
                for l in f:
                    l=l.strip()
                    l=list(map(conv,l.split(",")[:-1]))
                    data[m][algo].append(l)

    for m,mv in data.items():
        print("-----",m,"------")
        for algo,d in mv.items():
            alist=[0 for _ in range(2,35)]
            inf_loop=[0 for _ in range(2,35)]
            for l in d:
                alist[l[0]-2]+=l[3]
                if l[2]<29.9 and l[3]==0:
                    inf_loop[l[0]-2]+=1
            #print(algo,": ")
            print(alist)
            #print(inf_loop)
