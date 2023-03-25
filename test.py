import subprocess, os,glob, time
processPool=[]
exe = "./debug/CCBS"
map_address ="./Examples/test_case/case{}/map.graph"
task_address="./Examples/test_case/case{}/task.task"
output_address="./test_result.txt"

try:
    os.remove("test_result.txt")
    print("removed successfully")
except OSError as e:
    print("{e}")
for m in [1,2,"3a","3b","4a","4b",5,6,7,8,9,10,11,12]:
#for m in [2]:
    for ES in [False]:#,True]:
        for cr in [False,True]:
            for ds in [False,True]:
                cmd=[exe,"-m",map_address.format(m),
                     "-t",task_address.format(m),
                    "--HI_h","0","precision","0.00001",
                     "--debug","0",
                     "-o", output_address,
                     "-a",str(4.5),
                     "--agent_num",str(10),
                     "timelimit","30",
                     "--extra_info",str(m if isinstance(m,int) else int(m[0]))]
                if ES:
                    cmd+=["--ES"]
                if cr:
                    cmd+=["CR"]
                if ds:
                    cmd+=["--DS"]

                print(subprocess.list2cmdline(cmd))
                if (len(processPool)>=8):
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
