#include "simpleparser.h"
#include <fstream>
#include <sstream>
#include <iostream>

SimpleParser::SimpleParser(const std::string &fileName)
{
    std::fstream f(fileName.c_str(), std::fstream::in);
    if(!f.fail())
    {
        while(!f.eof())
        {
            char str[256];
            f.getline(str, 256);
            std::string line(str);

            if(line.empty()) continue;
            if(line[0] == '#') continue;

            size_t dot2pos = line.find_first_of(':');
            if(dot2pos != std::string::npos)
            {
                std::string param, values;
                param = line.substr(0, dot2pos);
                values = line.substr(dot2pos+1);

                std::stringstream ss(values);
                while(!ss.eof())
                {
                    double value;
                    ss >> value;
                    _task[param].push_back(value);
                }
            }
        }
        f.close();

        /*
        for(auto &item : _task)
        {
            std::cout << item.first << " {";
            for(auto &value : item.second) std::cout << value << ", ";
            std::cout << "}\n";
        }
        */
    }
}

void SimpleParser::GetTask(taskType &task)
{
    task = _task;
}
