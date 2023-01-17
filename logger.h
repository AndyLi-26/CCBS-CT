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

  public:
    logger(Config cfg) :config(cfg) {};
    void write_to_log_summary(const Solution &solution);
    void write_to_log_path(const Solution &solution, const Map &map);
    void write_nodes(const Map &map);
};

#endif // LOGGER_H
