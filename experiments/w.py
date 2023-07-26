import json
x={
        "r":"0.353553391",#"0.5 1.5 2.5 3.5 4.5",
        "m":"sparse dense super-dense",
        "es":"0",
        "ct":"0 1 2",
        "icp":"0 1",
        "ds":"1 0",
        "a":" ".join(map(str,range(2,41))),
        "i":" ".join(map(str,range(1,26)))
        }
y=json.dumps(x,indent=1)
print(y)
with open("./config.json","w") as f:
    print(y,file=f)
