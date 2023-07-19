import subprocess, os,glob, time,json
processPool=[]
exe = "../debug/CCBS"
map_address ="../Instances/roadmaps/{}/map.graph"
task_address="../Instances/roadmaps/{}/large_agent/{}_task.task"
output_address="{}-{}-{}-{}-{}.csv"
with open("./config.json","r") as f:
    config=json.loads(f.read())
    for k,v in config.items():
        config[k]=v.split(" ")
    print(config)

for r in config['r']:
    for m in config['m']:
        for es in config['es']:
            for ct in config['ct']:
                for ds in config['ds']:
                    for icp in config['icp']:
                        for a in config['a']:
                            for i in config['i']:
                                es_tag="es" if es=='1' else '0'
                                if ct=='0': ct_tag="0" 
                                elif ct=='1': ct_tag="ct"
                                elif ct=="2": ct_tag="ct_abs"
                                ds_tag="ds" if ds=='1' else '0'
                                icp_tag="icp" if icp=='1' else '0'

                                cmd=[exe,"-m",map_address.format(m),
                                     "-t",task_address.format(m,i),
                                    "--HI_h","2",
                                     "-o", output_address.format(m,es_tag,ct_tag,ds_tag,icp_tag),
                                     "-a",r,
                                     "--agent_num",a,
                                     "--timelimit","30",
                                     "--extra_info",i
                                ]
                                if es=='1':
                                    cmd+=["--ES"]
                                if ct=='1':
                                    cmd+=["--CT"]
                                elif ct=='2':
                                    cmd+=["--CT_abs"]
                                if ds=='1':
                                    cmd+=["--DS"]
                                if icp=='1':
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
                                
