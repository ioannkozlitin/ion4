#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <vector>

void outputArray(std::ostream& os, const std::string dataName, const std::vector<double> &data);
void outputTable(std::ostream& os, std::string tableName, const std::vector<std::vector<double>> &table);

#endif // OUTPUT_H
