import subprocess, os,glob, time
processPool=[]
exe = "./debug/CCBS"
map_address ="./Examples/test_case/case{}/map.graph"
task_address="./Examples/test_case/case{}/task.task"
output_address="./test_result.txt"
if os.path.exists(output_address):
    os.remove(output_address)
    print("File deleted successfully.")
else:
    print("File does not exist.")
for m in [1,2,"3a","3b","4a","4b",5,6,7,8,9,10,11,12]:
    #for m in [7,9,10,11]:
#for m in [2]:
    for ES in [False]:
        for ct in [0,1,2]:
            for ds in [False,True]:
                for icp in [False,True]:
                    cmd=[exe,"-m",map_address.format(m),
                         "-t",task_address.format(m),
                        "--HI_h","0","precision","0.00001",
                         "-o", output_address,
                         "-a",str(0.5),
                         "--agent_num",str(10),
                         "timelimit","30",
                         "--extra_info",str(m if isinstance(m,int) else int((m[0] + str(["a","b"].index(m[1])))))]
                    if ES:
                        cmd+=["--ES"]
                    if ct==1:
                        cmd+=["--CT"]
                    elif ct==2:
                        cmd+=["--CT_abs"]
                    if ds:
                        cmd+=["--DS"]
                    if icp:
                        cmd+=["--ICP"]

                    print(subprocess.list2cmdline(cmd))
                    if (len(processPool)>=4):
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
