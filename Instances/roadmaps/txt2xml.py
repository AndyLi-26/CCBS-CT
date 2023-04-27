m="./sparse/24task.task"
template='''    <agent start_id="{}" goal_id="{}"/>'''
head='''<?xml version="1.0" ?>
<root>'''
foot="</root>"
t=[]
with open(m,'r') as f:
    n=int(f.readline().strip())
    for i in range(n):
        n1,n2=map(int,f.readline().strip().split(' '))
        t.append((n1,n2))

with open("my_task.xml",'w') as f:
    print(head,file=f)
    for i in t:
        print(template.format(*i),file=f)
    print(foot,file=f)

