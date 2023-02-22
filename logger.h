#ifndef LOGGER_H
#define LOGGER_H
#include "const.h"
#include "structs.h"
#include "map.h"
#include <iostream>
#include <string>
#include <fstream>
class logger
{
  private:
    Config config;
    Solution solution;
  public:
    logger(Config cfg, Solution sol) :config(cfg),solution(sol) {};
    void write_to_log_summary();
    void write_to_log_path(const Map &map);
    void write_nodes(const Map &map);
    void write_exp_result(int task_ind);
};

#endif // LOGGER_H
