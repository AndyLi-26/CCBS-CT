import subprocess, os,glob, time,json
processPool=[]
exe = "../release/CCBS"
map_address ="../Instances/roadmaps/{}/map.graph"
task_address="../Instances/roadmaps/{}/large_agent/{}_task.task"
output_address="{}-{}-{}-{}.csv"
with open("./config.json","r") as f:
    config=json.loads(f.read())
    for k,v in config.items():
        config[k]=v.split(" ")
    print(config)
algo=[]
for r in config['r']:
    for m in config['m']:
        for es in config['es']:
            for ct in config['ct']:
                for ds in config['ds']:
                    for i in config['i']:
                        for a in config['a']:
                            es_tag="es" if es=='1' else '0'
                            if ct=='0': ct_tag="0"
                            elif ct=='1': ct_tag="ct"
                            elif ct=="2": ct_tag="ct_abs"
                            elif ct=="3": ct_tag="icp"
                            ds_tag="ds" if ds=='1' else '0'

                            cmd=[exe,"-m",map_address.format(m),
                                 "-t",task_address.format(m,i),
                                "--HI_h","2",
                                 "-o", output_address.format(m,es_tag,ct_tag,ds_tag),
                                 "-a",r,
                                 "--agent_num",a,
                                 "--timelimit","30",
                                 "--extra_info",i,
                                 "--Cardinal"
                            ]
                            if es=='1':
                                cmd+=["--ES"]
                            if ct=='1':
                                cmd+=["--CT"]
                            elif ct=='2':
                                cmd+=["--CT_abs"]
                            elif ct=='3':
                                cmd+=["--CT_abs"]
                                cmd+=["--ICP"]
                            if ds=='1':
                                cmd+=["--DS"]
                            print(subprocess.list2cmdline(cmd))
                            assert False

                            if (len(processPool)>=3):
                                finish = False
                                while not finish:
                                    time.sleep(1)
                                    for p in range(0,len(processPool)):
                                        if p >= len(processPool):
                                            break
                                        result=processPool[p].poll()
                                        if result is not None:
                                            if result!=0:
                                                print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                                                print(processPool[p].args)
                                                print(processPool[p].stderr)
                                                with open("errors.txt",'a') as f:
                                                    print(' '.join(processPool[p].args),file=f)
                                            processPool.pop(p)
                                            finish = True
                                            p-=1
                            else:
                                for p in range(0,len(processPool)):
                                        if p >= len(processPool):
                                            break

                                        result=processPool[p].poll()
                                        if result is not None:
                                            print("result",result)
                                            if result!=0:
                                                print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                                                print(processPool[p])
                                                print(processPool[p].args)
                                                with open("errors.txt",'a') as f:
                                                    print(' '.join(processPool[p].args),file=f)
                                            processPool.pop(p)
                                            finish = True
                                            p-=1

                            try:
                                processPool.append(subprocess.Popen(cmd))
                            except:
                                print(len(processPool))

