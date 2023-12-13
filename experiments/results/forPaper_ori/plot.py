import matplotlib.pyplot as plt
import os,sys
import numpy as np
global TOTAL
TOTAL=0
def align(d):
    d = sorted(d, key=lambda x: x[0])
    d = sorted(d, key=lambda x: x[6])
    d = sorted(d, key=lambda x: x[5])
    return d

def conv(s):
    try:
        return int(s)
    except:
        try:
            return float(s)
        except:
            return s

def readF(f):
    DATA=[]
    for l in f:
        temp=l.strip().split(',')
        temp=temp[:-1]
        temp=list(map(conv,temp))
        #print(temp)
        DATA.append(temp)
    return DATA

def getSuccRate(d):
    temp=set()
    for l in d:
        temp.add(l[1])
    x=list(temp)
    x.sort()
    y=[0 for _ in x]
    temp=[0 for _ in x]

    for l in d:
        idx=x.index(l[1])
        y[idx]+=l[2]
        temp[idx]+=1

    assert all(i == 25 for i in temp)
    #y = [i/25*100 for i in y]
    return x,y

def getSolvedInsT(x,d):
    y=[]
    for i in x:
        num=0
        for l in d:
            if l[7]<i: num+=l[8]
        y.append(num)
    return y

def getSolvedInsE(x,d):
    y=[]
    for i in x:
        num=0
        for l in d:
            if l[12]<i: num+=l[8]
        y.append(num)
    return y

def plotDATA(d,t,n):
# Create a figure and axis
    fig, ax = plt.subplots()
    for fn in d:
        ax.plot(d[fn][0], d[fn][1], label=fn)

# Add labels and a legend
    ax.set_xlabel('Agents')
    ax.set_ylabel('success rate')
    ax.set_title(t)
    ax.legend()
# Show the plot
    plt.show()

def comparsionPlot(D,t,n):
    idx=1
    for i in D:
        plt.subplot(3, 4, idx)
        idx+=1
        x1,y1=getSuccRate(D[i][0])
        x2,y2=getSuccRate(D[i][1])
        plt.plot(x1,y1,label=n[0])
        plt.plot(x2,y2,label=n[1])
        plt.title(i)
        plt.legend()

    plt.suptitle(t)
    plt.tight_layout()
    plt.show()

def countDead(d):
    c=0
    for l in d:
        c+=l[-1]
    return c


if __name__=="__main__":
    current_directory = os.getcwd()
    csv_files = [file for file in os.listdir(current_directory) if file.endswith('.csv')]
    DATAALL=[dict() for _ in 'asd']
    for csv_file in csv_files:
        file_path = os.path.join(current_directory, csv_file)
        with open(file_path, 'r') as f:
            #print(file_path)
            d=readF(f)
            if "sparse" in csv_file: i=0
            elif "super" in csv_file: i=2
            else: i=1
            k=csv_file.replace('.csv','')
            k=k.replace('sparse-','')
            k=k.replace('super-dense-','')
            k=k.replace('dense-','')

            DATAALL[i][k]=align(d)

    ms=["sparse",'dense',"super-dense"]
    solveIns=['0-0-0','0-0-ds','0-ct_abs-ds','0-icp-ds']
    fig, ax = plt.subplots()
    for i,d in enumerate(DATAALL):
        #print solved ins inT
        ys=[]
        plt.subplot(3, 1, i+1)
        for k in solveIns:
            print(k)
            x = np.arange(0, 30.01, 0.01).tolist()
            x= np.logspace(np.log10(0.0001), np.log10(30), num=100).tolist()
            y=getSolvedInsT(x,d[k])
            ys.append(y)
            plt.plot(y,x, label=k)

# Add labels and a legend
        plt.ylabel('time')
        plt.xlabel('solved Instances')
        plt.yscale('log')
        plt.title(ms[i])
        plt.legend()
# Show the plot
    plt.tight_layout()
    plt.show()

    plt.figure()
    for i,d in enumerate(DATAALL):
        #print solved ins in EXP
        ys=[]
        plt.subplot(3, 1, i+1)
        maxEXP=0
        for k in solveIns:
            for l in d[k]:
                maxEXP=max(maxEXP,l[12])

        print(maxEXP)
        for k in solveIns:
            print(k)
            #x = np.arange(0, 30.01, 0.01).tolist()
            x= np.logspace(np.log10(1), np.log10(maxEXP), num=100).tolist()
            y=getSolvedInsE(x,d[k])
            ys.append(y)
            plt.plot(y, x, label=k)

# Add labels and a legend
        plt.ylabel('Node expansion')
        plt.xlabel('solved Instances')
        plt.yscale('log')
        plt.title(ms[i])
        plt.legend()
# Show the plot
    plt.tight_layout()
    plt.show()











