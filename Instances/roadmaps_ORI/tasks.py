import xml.etree.ElementTree as ET
def xml2data(f):
    with open(f,'r') as f:
        s=f.read().strip()
    root = ET.fromstring(s)
    data = [(int(agent.attrib['start_id']), int(agent.attrib['goal_id'])) for agent in root.findall('agent')]
    return data
for m in ["sparse","dense","super-dense"]:
    for i in range(1,26):
        rf=f"./{m}/{i}_task.xml"
        d=xml2data(rf)
        
        wf=f"./{m}/{i}_task.task"
        with open(wf,"w") as f:
            print(len(d),file=f)
            [print(*l,file=f) for l in d]
