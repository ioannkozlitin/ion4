#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>

class Element
{
public:
    Element(unsigned int z);

    unsigned int Z; //Заряд элемента
    std::vector<double> fi; //Массив потенциалов ионизации в а.е.
};

double GetA(unsigned int z);

#endif // ELEMENT_H
