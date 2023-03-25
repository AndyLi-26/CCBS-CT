import sys
def remove_rubbish(s):
    x=s.index("<node")
    s=s[x:]
    s=s.strip()
    return s
    
def parserNode(info):
    b=info.index('>')
    b1=info[:b+1]
    b2=info[b+1:].strip()
    i=int(b1[11:-2])
    coord=tuple(map(float,b2[17:-7].split(',')))
    return coord

def parserEdge(info):
    b=info.index('>')
    b1=info[:b+1]
    b2=info[b+1:].strip()
    b1=b1.split('"n')[1:]
    parse=lambda s:int(s[:s.index('"')])
    b1=tuple(map(parse,b1))
    return b1
    
if __name__=="__main__":
    for name in ["10","11","12"]:
        with open("case"+name+"/map.xml",'r') as f:
            s=f.read().strip()
        s=remove_rubbish(s)
        nodes=list(map(lambda x:x.strip(),s.split("</node>")))
        edges=list(map(lambda x:x.strip(),nodes[-1].split("</edge>")))
        nodes=nodes[:-1]
        edges=edges[:-1]
        print(nodes[0])
        print(edges[0])
        newName=("case"+name+"/map.xml").split(".")
        newName=newName[:-1]+["graph"]
        newName=".".join(newName)
        with open(newName,'w') as f:
            print(len(nodes),file=f)
            for info in nodes:
                print(*parserNode(info),file=f)
            print(len(edges),file=f)
            for info in edges:
                print(*parserEdge(info),file=f)
