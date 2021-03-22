#ifndef SIMPLEPARSER_H
#define SIMPLEPARSER_H

#include <map>
#include <vector>
#include <string>

typedef std::map<std::string, std::vector<double>> taskType;

class SimpleParser
{
public:
    SimpleParser(const std::string &fileName);
    void GetTask(taskType &task);

private:
    taskType _task;
};

#endif // SIMPLEPARSER_H
