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
        temp.add(l[0])
    x=list(temp)
    x.sort()
    y=[0 for _ in x]
    temp=[0 for _ in x]

    for l in d:
        idx=x.index(l[0])
        y[idx]+=l[8]
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

def plotTvSolved(ms,solveIns,markers):
    #plt.figure()
    #fig, ax = plt.subplots()
    for i,d in enumerate(DATAALL):
        #print solved ins inT
        ys=[]
        plt.subplot(3, 3, i*3+1)
        for j,k in enumerate(solveIns):
            print(k)
            x = np.arange(0, 30.01, 0.01).tolist()
            x= np.logspace(np.log10(0.0001), np.log10(30), num=30).tolist()
            y=getSolvedInsT(x,d[k])
            ys.append(y)
            plt.plot(y,x, label=name[j],marker=markers[j],markerfacecolor="None",alpha=0.5,dashes=ls[j])
        if i==0: plt.legend(loc='lower center', bbox_to_anchor=(1.6, -3.1), ncol=7)

# Add labels and a legend
        plt.ylabel('time')
        plt.xlabel('solved Instances')
        plt.yscale('log')
        plt.title(ms[i])
# Show the plot
    #plt.tight_layout()

def plotEXPvSolved(ms,solveIns,markers):
    #plt.figure()
    for i,d in enumerate(DATAALL):
        #print solved ins in EXP
        ys=[]
        plt.subplot(3, 3, i*3+2)
        maxEXP=0
        for k in solveIns:
            for l in d[k]:
                maxEXP=max(maxEXP,l[12])

        print(maxEXP)
        for j,k in enumerate(solveIns):
            print(k)
            #x = np.arange(0, 30.01, 0.01).tolist()
            x= np.logspace(np.log10(1), np.log10(maxEXP), num=30).tolist()
            y=getSolvedInsE(x,d[k])
            ys.append(y)
            plt.plot(y, x, label=name[j],marker=markers[j],markerfacecolor="None",alpha=0.5,dashes=ls[j])

# Add labels and a legend
        plt.ylabel('Node expansion')
        plt.xlabel('solved Instances')
        plt.yscale('log')
        plt.title(ms[i])
# Show the plot
    #plt.tight_layout()

def plotSucc(ms,solveIns,markers):
    #plt.figure()
    for i,d in enumerate(DATAALL):
        #print solved ins in EXP
        plt.subplot(3, 3, i*3+3)

        for j,k in enumerate(solveIns):
            print(k)
            #x = np.arange(0, 30.01, 0.01).tolist()
            #x= np.logspace(np.log10(1), np.log10(maxEXP), num=100).tolist()
            x,y=getSuccRate(d[k])
            plt.plot(x, y, label=name[j],marker=markers[j],markerfacecolor="None",alpha=0.5,dashes=ls[j])

# Add labels and a legend
        plt.ylabel('Agents')
        plt.xlabel('solved Instances')
        plt.title(ms[i])
# Show the plot
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.15,hspace=0.35)
    plt.show()

def plotTvT():
    for i,c in enumerate(comp):
        plt.subplot(2, 2, i+3)
        x=[]
        y=[]
        solved=[0,0,0,0]#we solved, they solved, both solve, no solved
        for i,d in enumerate(DATAALL):
            algo1=solveIns[c[0]]
            n1=name[c[0]]
            d1=d[algo1]
            algo2=solveIns[c[1]]
            n2=name[c[1]]
            d2=d[algo2]
            for j in range(len(d1)):
                solved1=d1[j][8]
                solved2=d2[j][8]

                if not solved1 and solved2:
                    solved[0]+=1
                elif solved1 and not solved2:
                    solved[1]+=1
                elif solved1 and solved2:
                    solved[2]+=1
                else:
                    solved[3]+=1

                if solved1:
                    x.append(d1[j][7])
                else:
                    x.append(30)

                if solved2:
                    y.append(d2[j][7])
                else:
                    y.append(30)
        c=[]
        result=[0,0,0]#better, same, worse
        for i in range(len(x)):
            if x[i]==30 and y[i]==30:
                c.append('black')
            elif (x[i]==30 and y[i]!=30) or x[i]>y[i]+0.5:
                c.append('green')
                result[0]+=1
            elif (y[i]==30 and x[i]!=30) or y[i]>x[i]+0.5:
                c.append('red')
                result[2]+=1
            else:
                c.append('blue')
                result[1]+=1


        ref= np.logspace(np.log10(0.00005), np.log10(30), num=110).tolist()
        ref100x= np.logspace(np.log10(0.005), np.log10(30), num=110).tolist()
        ref100y= np.logspace(np.log10(0.00005), np.log10(0.3), num=110).tolist()
        ref10x= np.logspace(np.log10(0.0005), np.log10(30), num=110).tolist()
        ref10y= np.logspace(np.log10(0.00005), np.log10(3), num=110).tolist()
        ref2x= np.logspace(np.log10(0.0001), np.log10(30), num=110).tolist()
        ref2y= np.logspace(np.log10(0.00005), np.log10(15), num=110).tolist()
        plt.plot(ref,ref,c='blue')
        plt.plot(ref2x,ref2y,c="blue",dashes=(10,10))
        plt.text(ref2x[-1], ref2y[-1]-1, "2x", color='blue', fontsize=15, ha='right', va='top',fontweight="bold",rotation=40,rotation_mode="default",snap=True)
        plt.plot(ref10x,ref10y,c="blue",dashes=(10,10))
        plt.text(ref10x[-1], ref10y[-1]-0.2, "10x", color='blue', fontsize=15, ha='right', va='top',fontweight="bold",rotation=35,rotation_mode="default",snap=True)
        plt.plot(ref100x,ref100y,c="blue",dashes=(10,10))
        plt.text(ref100x[-1], ref100y[-1]-0.05, "100x", color='blue', fontsize=15, ha='right', va='top',fontweight="bold",rotation=30,rotation_mode="default",snap=True)
        plt.scatter(x,y,c=c)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel(n1)
        plt.ylabel(n2)
        plt.title("runtime VS runtime \n{} better, {} same, {} worse, we solved {} more instances".format(*result,solved[0]-solved[1]))

def plotExpvExp():
    for i,c in enumerate(comp):
        plt.subplot(2, 2, i+1)
        m1=m2=0

        for i,d in enumerate(DATAALL):
            algo1=solveIns[c[0]]
            d1=d[algo1]
            algo2=solveIns[c[1]]
            d2=d[algo2]
            m1=max(max(d1,key=lambda x:x[12])[12],m1)
            m2=max(max(d2,key=lambda x:x[12])[12],m2)
        print(m1,m2)
        m1=int(m1*1.1)
        m2=int(m2*1.1)

        x=[]
        y=[]
        solved=[0,0,0,0]#we solved, they solved, both solve, no solved
        for i,d in enumerate(DATAALL):
            algo1=solveIns[c[0]]
            n1=name[c[0]]
            d1=d[algo1]
            algo2=solveIns[c[1]]
            n2=name[c[1]]
            d2=d[algo2]
            for j in range(len(d1)):
                solved1=d1[j][8]
                solved2=d2[j][8]
                if not solved1 and solved2:
                    solved[0]+=1
                elif solved1 and not solved2:
                    solved[1]+=1
                elif solved1 and solved2:
                    solved[2]+=1
                else:
                    solved[3]+=1

                if d1[j][8]==1:
                    x.append(d1[j][12])
                else:
                    x.append(m1)

                if d2[j][8]==1:
                    y.append(d2[j][12])
                else:
                    y.append(m2)
        c=[]
        result=[0,0,0]#better, same, worse
        for i in range(len(x)):
            if x[i]==m1 and y[i]==m2:
                c.append('black')
            elif x[i]==y[i]:
                c.append('blue')
                result[1]+=1
            elif x[i]>y[i]:
                c.append('green')
                result[0]+=1
            else:
                c.append('red')
                result[2]+=1
        ref= np.logspace(np.log10(1), np.log10(max(m1,m2)), num=110).tolist()
        ref2x= np.logspace(np.log10(2), np.log10(max(m1,m2)), num=110).tolist()
        ref2y= np.logspace(np.log10(1), np.log10(max(m1,m2)/2), num=110).tolist()
        ref10x= np.logspace(np.log10(10), np.log10(max(m1,m2)), num=110).tolist()
        ref10y= np.logspace(np.log10(1), np.log10(max(m1,m2)/10), num=110).tolist()
        ref100x= np.logspace(np.log10(100), np.log10(max(m1,m2)), num=110).tolist()
        ref100y= np.logspace(np.log10(1), np.log10(max(m1,m2)/100), num=110).tolist()
        plt.plot(ref,ref,c='blue')
        plt.plot(ref2x,ref2y,c="blue",dashes=(10,10))
        plt.text(ref2x[-1], ref2y[-1]-20000, "2x", color='blue', fontsize=15, ha='right', va='top',fontweight="bold",rotation=40,rotation_mode="default",snap=True)
        plt.plot(ref10x,ref10y,c="blue",dashes=(10,10))
        plt.text(ref10x[-1], ref10y[-1]-5000, "10x", color='blue', fontsize=15, ha='right', va='top',fontweight="bold",rotation=35,rotation_mode="default",snap=True)
        plt.plot(ref100x,ref100y,c="blue",dashes=(10,10))
        plt.text(ref100x[-1], ref100y[-1]-500, "100x", color='blue', fontsize=15, ha='right', va='top',fontweight="bold",rotation=30,rotation_mode="default",snap=True)
        plt.scatter(x,y,c=c)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel(n1)
        plt.ylabel(n2)
        plt.title("Expansion VS Expansion \n {} better, {} same, {} worse, we solved {} more instances".format(*result,solved[0]-solved[1]))



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
    solveIns=['0-0-ds-0','0-0-ds-2','0-ct_abs-ds-0','0-ct_abs-ds-2','0-icp-ds-0','0-icp-ds-2']
    name=['DS','DSHP','DS-CT','DSHP-CT','DS-GCT','DSHP-GCT']
    markers=['.','o','v','^','1','s','*']
    ls=[(),(1, 1),(1, 5), (5, 8), (3, 10, 1, 10), (3, 5, 1, 5), (3, 5, 2, 5, 1, 5)]
    plotTvSolved(ms,solveIns,markers)
    plotEXPvSolved(ms,solveIns,markers)
    plotSucc(ms,solveIns,markers)
    comp=[(1,5),(0,4)]
    plt.figure()
    plotExpvExp()
    plotTvT()
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.15,hspace=0.3)
    plt.show()
