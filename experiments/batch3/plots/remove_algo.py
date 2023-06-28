import os
current_folder = os.getcwd()
files = os.listdir(current_folder)
# Iterate over each file and print its name if it's a CSV file
file_list=[]
for file_name in files:
    if file_name.endswith('.csv'):
        file_list.append(file_name)

for fn in file_list:
    data=[]
    with open(fn) as f:
        for l in f:
            l=l.strip().split(',')
            l=[l[0]]+l[5:]
            l=",".join(l)
            data.append(l)

    with open(fn,'w') as f:
        [print(i,file=f) for i in data]

