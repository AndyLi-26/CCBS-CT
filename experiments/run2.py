import subprocess, os,glob, time,json
processPool=[]
exe = "../debug/CCBS"
map_address ="../Instances/roadmaps/{}/map.graph"
task_address="../Instances/roadmaps/{}/large_agent/{}_task.task"
output_address="{}-{}-{}-{}-{}.csv"
algo=["","ds","ds ct1","ds ct2","ds icp","ds icp ct2"]
for m in ["sparse","dense","super-dense"]:
    for a in range(3,41):
        for i in range(1,26):
            for alg in algo:

                es_tag="es" if "es" in alg else '0'
                if "ct1" in alg: ct_tag="ct"
                elif "ct2" in alg: ct_tag="ct_abs"
                else: ct_tag="0"
                ds_tag="ds" if "ds" in alg else '0'
                icp_tag="icp" if "icp" in alg else '0'

                cmd=[exe,"-m",map_address.format(m),
                     "-t",task_address.format(m,i),
                    "--HI_h","2",
                     "-o", output_address.format(m,es_tag,ct_tag,ds_tag,icp_tag),
                     "-a","0.353553391",
                     "--agent_num",str(a),
                     "--timelimit","30",
                     "--extra_info",str(i)
                ]

                if "es" in alg:
                    cmd+=["--ES"]
                if "ct1" in alg:
                    cmd+=["--CT"]
                elif "ct2" in alg:
                    cmd+=["--CT_abs"]
                if "ds" in alg:
                    cmd+=["--DS"]
                if "icp" in alg:
                    cmd+=["--ICP"]
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
                
                try:
                    processPool.append(subprocess.Popen(cmd))
                except:
                    print(len(processPool))
                
