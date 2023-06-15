import json
x={
        "r":"0.35355",#"0.5 1.5 2.5 3.5 4.5",
        "m":"sparse dense super-dense",
        "es":"0",
        "cr":"0",
        "card": "1 0",
        "ds":"1 0",
        "a":" ".join(map(str,range(2,35))),
        "i":" ".join(map(str,range(1,26)))
        }
y=json.dumps(x,indent=1)
print(y)
with open("./config.json","w") as f:
    print(y,file=f)
