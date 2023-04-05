import subprocess, os,glob, time
processPool=[]
exe = "../debug/CCBS"
map_address ="../Instances/roadmaps/{}/map.graph"
task_address="../Instances/roadmaps/{}/{}task.task"
output_address="./{}.csv"

for a_size in [4.5]:#[0.5,1.5,2.5,3.5,4.5]:
    for m in ["sparse"]:#["sparse","dense","super-dense"]:
        for e_sp in [False]:#[False,True]:
            for cr in [False,True]:
                for ds in [False,True]:
                    for a in range(5,6):#(5,26):
                        for i in range(1,26):
                            cmd=[exe,"-m",map_address.format(m),
                                 "-t",task_address.format(m,i),
                                "--HI_h","0","precision","0.00001",
                                 "-o", output_address.format(m),
                                 "-a",str(a_size),
                                 "--agent_num",str(a),
                                 "timelimit","30",
                                 "--extra_info",str(i)
                            ]
                            if e_sp:
                                cmd+=["--ES"]
                            if cr:
                                cmd+=["--CR"]
                            if ds:
                                cmd+=["--DS"]
                            print(subprocess.list2cmdline(cmd))
                            if (len(processPool)>=9):
                                finish = False
                                while not finish:
                                    time.sleep(1)
                                    for p in range(0,len(processPool)):
                                        if p >= len(processPool):
                                            break
                                        if processPool[p].poll() is not None:
                                            processPool.pop(p)
                                            finish = True
                                            p-=1
                            else:
                                for p in range(0,len(processPool)):
                                        if p >= len(processPool):
                                            break
                                        if processPool[p].poll() is not None:
                                            processPool.pop(p)
                                            finish = True
                                            p-=1
                            
                            print(subprocess.list2cmdline(cmd))
                            try:
                                processPool.append(subprocess.Popen(cmd))
                            except:
                                print(len(processPool))
                            print(processPool)
