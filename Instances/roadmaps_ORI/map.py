import xml.etree.ElementTree as ET
with open("./super-dense/map.xml",'r') as f:
    s=f.read().strip()
root = ET.fromstring(s)[2]
nodes = []
pref="{http://graphml.graphdrawing.org/xmlns}"
for n in root.findall(pref+"node"):
    data_string = n.find(pref+'data').text
    numbers = tuple(float(num) for num in data_string.split(','))
    nodes.append(numbers)
edges = [(int(e.attrib['source'].strip('n')), int(e.attrib['target'].strip('n'))) for e in root.findall(pref+'edge')]
with open("./super-dense/map.graph","w") as f:
    print(len(nodes),file=f)
    [print(*i,file=f) for i in nodes]
    print(len(edges),file=f)
    [print(*i,file=f) for i in edges]
