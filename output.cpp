#include "output.h"

void outputArray(std::ostream &os, const std::string dataName, const std::vector<double> &data)
{
    os << dataName << " = [";
    for(std::vector<double>::const_iterator it = data.begin(); it!= data.end(); ++it) {os << *it << " ";};
    os << "];" << std::endl;
}

void outputTable(std::ostream &os, std::string tableName, const std::vector<std::vector<double> > &table)
{
    os << tableName << " = [" << std::endl;
    for (size_t i = 0; i < table.size(); ++i)
    {
        const std::vector<double>& line = table[i];
        for (size_t j = 0; j < line.size(); ++j)
        {
            os << (line[j]) << " ";
        }
        os << std::endl;
    }
    os << "];" << std::endl;
}
