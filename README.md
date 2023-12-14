# CCBS-CT

clone the repo
```
git clone https://github.com/PathPlanning/Continuous-CBS.git
```
or direct downloading.

```bash
mkdir debug
cd debug
cmake ../
make -j
```
## Options
Here are some CLI parameter you can use, full list is in main.cpp, there is also more parameter can be modified in cons.h:
* `-m <map file location>` - the map file
* `-t <task file location>` - the task file
* `agent_num <number of agent>` - number of agent in int
* `-a <agent size in float>` -agent size in float
* `timelimit` - timelimit

## Launch
Here is an example of minimum information required to launch the code
```
../debug/CCBS -m ../Instances/roadmaps/sparse/map.graph -t ../Instances/roadmaps/sparse/ori_set/10_task.task --agent_num 3
```
Here is an example of more advance constrol on the code
```
../debug/CCBS -m ../Instances/roadmaps/sparse/map.graph -t ../Instances/roadmaps/sparse/ori_set/10_task.task --HI_h 0 -o test.csv -a 0.5 --agent_num 3 --timelimit 30 --extra_info 1 --DS --CT_abs --ICP
```
