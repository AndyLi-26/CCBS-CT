for i in range(1,26):
    tasks=[]
    print(i)
    with open(str(i)+"task.task",'r') as f:
        n=int(f.readline())
        for _ in range(n):
            t=f.readline().strip().split(' ')
            t=list(map(int,t))
            if t[0]==85: t[0]=120
            if t[1]==85: t[1]=120
            tasks.append(t)

    with open(str(i)+"task.task",'w') as f:
        print(len(tasks),file=f)
        [print(*i,file=f) for i in tasks ]

