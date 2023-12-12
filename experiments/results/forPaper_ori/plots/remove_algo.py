import os
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

current_folder = os.getcwd()
files = os.listdir(current_folder)
# Iterate over each file and print its name if it's a CSV file
file_list=[]
for file_name in files:
    if file_name.endswith('.csv'):
        file_list.append(file_name)

for fn in file_list:
    temp=[]
    with open(fn) as f:
        for l in f:
            l=l.strip().split(',')
            l=list(map(convert,l))
            if l[0]==2: continue
            temp.append(l)
        temp=s_l(temp)
        data=[] 
        for l in temp:
            l=[l[0]]+l[5:]
            l=",".join(map(str,l))
            data.append(l)

    with open(fn,'w') as f:
        [print(i,file=f) for i in data]

