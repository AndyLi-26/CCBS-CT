# CCBS-CT

clone the repo
```
git clone https://github.com/PathPlanning/Continuous-CBS.git
```
or direct downloading.

Use CMake and make to compile the code, you can create release or debug version with CMake tag
```bash
mkdir release
cd release
cmake ../
make -j
```
## Options
Here are some basic CLI parameter you can use, full list with description is in main.cpp, there is also more parameter can be modified in cons.h:
* `-m <map file location>` - the map file
* `-t <task file location>` - the task file
* `agent_num <number of agent>` - number of agent in int
* `-a <agent size in float>` - agent size in float (optional)
* `timelimit` - timelimit (optional)

## Launch
Here is an example of minimum information required to launch the code
```
./CCBS -m ../Instances/roadmaps/sparse/map.graph -t ../Instances/roadmaps/sparse/ori_set/10_task.task --agent_num 3
```
Here is an example of more advance control on the code
```
./CCBS -m ../Instances/roadmaps/sparse/map.graph -t ../Instances/roadmaps/sparse/ori_set/10_task.task --HI_h 0 -o test.csv -a 0.5 --agent_num 3 --timelimit 30 --extra_info 1 --DS --CT_abs --ICP
```
## Experiment
The script to run experiment is in the experiment folder
