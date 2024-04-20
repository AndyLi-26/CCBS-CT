import sys
def readF(f):
    nnum=int(f.readline().strip())
    for _ in range(nnum):
        l=f.readline().strip().split(" ")
        l=(float(l[0]),float(l[1]))
        N.append(l)

    enum=int(f.readline().strip())
    for _ in range(enum):
        l=f.readline().strip().split(" ")
        print(l)
        l=(int(l[0]),int(l[1]))
        E.append(l)

def writeF(f):
    sn='''    <node id="n{}">
      <data key="key0">{},{}</data>
    </node>'''
    se='''    <edge id="e{}" source="n{}" target="n{}">
      <data key="key1">1</data>
    </edge>'''
    beg='''<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <key id="key0" for="node" attr.name="coords" attr.type="string" />
  <key id="key1" for="edge" attr.name="weight" attr.type="double" />
  <graph id="G" edgedefault="directed" parse.nodeids="free" parse.edgeids="canonical" parse.order="nodesfirst">'''
    end='''  </graph>
</graphml>'''

    print(beg,file=f)
    for i,n in enumerate(N):
        print(sn.format(i,*n),file=f)

    for i,e in enumerate(E):
        print(se.format(i,*e),file=f)

    print(end,file=f)


if __name__ == "__main__":
    E=[]
    N=[]
    fn=sys.argv[1]
    with open(fn) as f:
        readF(f)
    print(N,E)



    outfn=".".join(fn.split(".")[:-1])+".xml"
    print(outfn)
    with open(outfn,"w") as f:
        writeF(f)







