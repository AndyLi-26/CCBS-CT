import subprocess, os,glob, time,json
processPool=[]
exe = "../debug/CCBS"
map_address ="../Instances/roadmaps_ORI/{}/map.graph"
task_address="../Instances/roadmaps_ORI/{}/{}_task.task"
output_address="{}-{}-{}-{}.csv"
with open("./config.json","r") as f:
    config=json.loads(f.read())
    for k,v in config.items():
        config[k]=v.split(" ")
    print(config)

for r in config['r']:
    for m in config['m']:
        for es in config['es']:
            for cr in config['cr']:
                for ds in config['ds']:
                    for card in config['card']:
                        for a in config['a']:
                            for i in config['i']:
                                es_tag="es" if es=='1' else '0'
                                cr_tag="cr" if cr=='1' else '0'
                                ds_tag="ds" if ds=='1' else '0'

                                cmd=[exe,"-m",map_address.format(m),
                                     "-t",task_address.format(m,i),
                                    "--HI_h","0",
                                     "-o", output_address.format(m,es_tag,cr_tag,ds_tag),
                                     "-a",r,
                                     "--agent_num",a,
                                     "--timelimit","120",
                                     "--extra_info",i
                                ]
                                if card=='1':
                                    cmd+=["--Cardinal"]
                                if es=='1':
                                    cmd+=["--ES"]
                                if cr=='1':
                                    cmd+=["--CR"]
                                if ds=='1':
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
                                
                                try:
                                    processPool.append(subprocess.Popen(cmd))
                                except:
                                    print(len(processPool))
                                
