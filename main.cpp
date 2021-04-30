#include "saha/src/elements.h"
#include "saha/src/sahasolver.h"
#include "saha/src/atom_ed.h"
#include "saha/saha.h"
#include <math.h>

#include "mix/mixdata.h"
#include "mix/sahamixsolver.h"
#include "raizer/raizermixsolver.h"

#include "saha/src/atom_ed.h"

#include "simpleparser.h"

#include <cstdio>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <map>

void outputArray(std::ostream& os, const std::string dataName, const std::vector<double> &data)
{
   os << dataName << " = [";
   for(std::vector<double>::const_iterator it = data.begin(); it!= data.end(); ++it) {os << *it << " ";};
   os << "];" << std::endl;
}

void outputTable(std::ostream& os, std::string tableName, const std::vector<std::vector<double>> &table)
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

void outputTable(std::ostream& os, const std::string &tableName, const std::vector<std::vector<SahaPoint>> &table, std::function<double(const SahaPoint&)> accessor, bool printT)
{
    os << tableName << " = [" << std::endl;
    for (size_t i = 0; i < table.size(); ++i)
    {
        const std::vector<SahaPoint>& line = table[i];
        if(printT)
        {
            os << log10(line[0].T * eFi) << " " << log10(line[0].t * eFi) << " ";
        }
        for (size_t j = 0; j < line.size(); ++j)
        {
            os << accessor(line[j]) << " ";
        }
        os << std::endl;
    }
    os << "];" << std::endl;
}

void CrashTest(double rCoeff, double lgVMin, double lgVMax, double lgVStep, double lgTMin, double lgTMax, double lgTStep)
{
	for (int Z = 1; Z <= 103; Z++)
	{
		const TElement elem(Z, rCoeff);
        SahaSolver solver(elem, 0);

		printf("[%d] ",Z);

		for (double lgT = lgTMax; lgT > lgTMin; lgT -= lgTStep)
		{
			for (double lgV = lgVMin; lgV < lgVMax; lgV += lgVStep)
			{
				SahaPoint res = solver.Calculate_lgTeV_lgVae(lgT, lgV);
			}
		}

	}
}

void calculator(unsigned int Z, double rCoeff, double lgVMin, double lgVMax, double lgVStep, double lgTMin, double lgTMax, double lgTStep, std::string filename)
{
    const TElement elem(Z, rCoeff, true);

    std::vector<double> lgTPhys;
    std::vector<double> lgVa;
    std::vector<double> lgRho;

    std::vector<std::vector<double>> ionizationTable;
    std::vector<std::vector<double>> pTable;
    std::vector<std::vector<double>> eTable;

    for (double lgT = lgTMax; lgT > lgTMin; lgT -= lgTStep) lgTPhys.push_back(lgT);

    double roConst = log10(eRo*elem.A);
    for (double lgV = lgVMin; lgV < lgVMax; lgV += lgVStep)
    {
        lgVa.push_back(lgV);
        lgRho.push_back(roConst - lgV);
    }

    for (double lgT = lgTMax; lgT > lgTMin; lgT -= lgTStep)
    {
        std::cout << "[" << lgT << "]" << std::flush;
        std::vector<double> ionizationLine;
        std::vector<double> pLine;
        std::vector<double> eLine;
        for (double lgV = lgVMin; lgV < lgVMax; lgV += lgVStep)
        {
            TElement elemCopy = elem;
            //elemCopy.softIon(pow(10,lgV));
            SahaSolver solver(elemCopy, 0);
            SahaPoint res = solver.Calculate_lgTeV_lgVae(lgT,lgV);
            ionizationLine.push_back(res.Xe);
            pLine.push_back(res.P);
            eLine.push_back(res.E);
        }
        ionizationTable.push_back(ionizationLine);
        pTable.push_back(pLine);
        eTable.push_back(eLine);
    }

    std::fstream f(filename.c_str(), std::fstream::out);
    f << std::scientific;

    f << "Z=" << Z << ";" << std::endl;
    outputArray(f, "lgT", lgTPhys);
    outputArray(f, "lgV", lgVa);
    outputArray(f, "lgRho", lgRho);
    outputTable(f, "xe_Saha", ionizationTable);
    outputTable(f, "P_Saha", pTable);
    outputTable(f, "E_Saha", eTable);
}

void calculatorRho_eV(unsigned int Z, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, std::string filename)
{
    double roConst = log10(eRo*elements::GetA(Z));
    double lgVMin = roConst - lgRhoMax;
    double lgVMax = roConst - lgRhoMin;
    double lgVStep = lgRhoStep;

    calculator(Z, rCoeff, lgVMin, lgVMax, lgVStep, lgTMin, lgTMax, lgTStep, filename);
}

void calculatorMix(const std::vector<unsigned int> &Z, const std::vector<double> &x, double rCoeff, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, double lgtDiffTmin, double lgtDiffTmax, std::string filename)
{
    SahaMixSolver mixSolver;

    std::vector<double> lgVa;
    std::vector<double> _lgRho;

    std::vector<std::vector<SahaPoint>> fullTable;
    std::vector<double> lgTPhys, TkeV;
    //std::vector<std::vector<std::vector<std::vector<double>>>> xxTable;

    //xxTable.resize(Z.size());
    //for(int i = 0; i < Z.size(); i++) xxTable[i].resize(Z[i] + 1);

    MixData md(Z, x, rCoeff, rCoeff > 0, rCoeff > 0, 1, 1, 1);
    bool is2T = (lgtDiffTmax > lgtDiffTmin);

    for (double lgT = lgTMax; lgT > lgTMin - lgTStep / 2.0; lgT -= lgTStep)
    for (double lgt = lgT + lgtDiffTmax; lgt > lgT + lgtDiffTmin - lgTStep / 2.0; lgt -= lgTStep)
    {        
        std::cout << "[" << lgt << " " << lgT << "]" << std::flush;
        std::vector<SahaPoint> fullLine;

        if(!is2T)
        {
            lgTPhys.push_back(lgT);
            TkeV.push_back(0.001 * pow(10, lgT));
        }

        //std::vector<std::vector<std::vector<double>>> xxLines;

        //xxLines.resize(Z.size());
        //for(int i = 0; i < Z.size(); i++) xxLines[i].resize(Z[i] + 1);

        bool fillFlag = _lgRho.empty();

        for (double lgRho = lgRhoMax; lgRho > lgRhoMin - lgRhoStep / 2.0; lgRho -= lgRhoStep)
        {
            //MixData md(Z, x, rCoeff, true, true, pow(10, lgT), pow(10, lgRho));
            md.SetTeVRho(pow(10, lgt), pow(10, lgT), pow(10, lgRho));

            mixSolver.GetFullIonizationInfo(md);
            fullLine.push_back(md.GetSahaPoint());

            //std::cout << "{" << xe - md.xe() << "}";

            /*for(int i = 0; i < Z.size(); i++)
            {
                for(int j = 0; j <= Z[i]; j++) xxLines[i][j].push_back(md.xx[i][j]);
            }*/

            if(fillFlag)
            {
                _lgRho.push_back(lgRho);
                lgVa.push_back(log10(md.GetFullV()));
            }
        }

        /*for(int i = 0; i < Z.size(); i++)
        {
            for(int j = 0; j <= Z[i]; j++)
            {
                xxTable[i][j].push_back(xxLines[i][j]);
            }
        }*/

        fullTable.push_back(fullLine);
    }

    std::fstream f(filename.c_str(), std::fstream::out);
    f << std::scientific;

    f << "Z=[";for(auto &z : Z) f << z << " ";f << "];" << std::endl;
    f << "x=[";for(auto &_x : x) f << _x << " ";f << "];" << std::endl;

    outputArray(f, "lgV", lgVa);
    outputArray(f, "lgRho", _lgRho);    
    if(!is2T)
    {
        outputArray(f, "lgT", lgTPhys);
        outputArray(f, "TkeV", TkeV);
    }

    outputTable(f, "xe", fullTable, std::mem_fn(&SahaPoint::Xe), is2T);
    outputTable(f, "mu", fullTable, std::mem_fn(&SahaPoint::M), is2T);
    outputTable(f, "F", fullTable, std::mem_fn(&SahaPoint::F), is2T);
    outputTable(f, "P", fullTable, std::mem_fn(&SahaPoint::P), is2T);
    outputTable(f, "Pi", fullTable, std::mem_fn(&SahaPoint::Pi), is2T);
    outputTable(f, "Pe", fullTable, std::mem_fn(&SahaPoint::Pe), is2T);
    outputTable(f, "E", fullTable, std::mem_fn(&SahaPoint::E), is2T);
    outputTable(f, "Ei", fullTable, std::mem_fn(&SahaPoint::Ei), is2T);
    outputTable(f, "Ee", fullTable, std::mem_fn(&SahaPoint::Ee), is2T);
    outputTable(f, "S", fullTable, std::mem_fn(&SahaPoint::S), is2T);
    outputTable(f, "Si", fullTable, std::mem_fn(&SahaPoint::Si), is2T);
    outputTable(f, "Se", fullTable, std::mem_fn(&SahaPoint::Se), is2T);

    /*for(int i = 0; i < Z.size(); i++)
    {
        for(int j = 0; j <= Z[i]; j++)
        {
            outputTable(f, "x_"+std::to_string(Z[i])+"_"+std::to_string(j), xxTable[i][j]);
        }
    }*/

}

void calculatorMixRaizer(const std::vector<unsigned int> &Z, const std::vector<double> &x, double lgRhoMin, double lgRhoMax, double lgRhoStep, double lgTMin, double lgTMax, double lgTStep, std::string filename)
{
    RaizerMixSolver mixSolver;

    std::vector<double> lgTPhys;
    std::vector<double> lgVa;
    std::vector<double> _lgRho;

    std::vector<std::vector<double>> ionizationTable;

    for (double lgT = lgTMax; lgT > lgTMin - lgTStep / 2.0; lgT -= lgTStep)
    {
        lgTPhys.push_back(lgT);
        std::cout << "[" << lgT << "]" << std::flush;
        std::vector<double> ionizationLine;

        bool fillFlag = _lgRho.empty();

        for (double lgRho = lgRhoMax; lgRho > lgRhoMin - lgRhoStep / 2.0; lgRho -= lgRhoStep)
        {
            MixData md(Z, x, 0.0, false, false, pow(10, lgT), pow(10, lgT), pow(10, lgRho));
            ionizationLine.push_back(mixSolver(md));

            if(fillFlag)
            {
                _lgRho.push_back(lgRho);
                lgVa.push_back(log10(md.GetFullV()));
            }
        }

        ionizationTable.push_back(ionizationLine);
    }

    std::fstream f(filename.c_str(), std::fstream::out);
    f << std::scientific;

    f << "Z=[";for(auto &z : Z) f << z << " ";f << "];" << std::endl;
    f << "x=[";for(auto &_x : x) f << _x << " ";f << "];" << std::endl;

    outputArray(f, "lgT", lgTPhys);
    outputArray(f, "lgV", lgVa);
    outputArray(f, "lgRho", _lgRho);
    outputTable(f, "xe_Saha", ionizationTable);
}

int main(int argc, char **argv)
{	
	try
	{
        std::string taskFileName, outputFileName;
        if(argc != 3)
        {
            std::cout << "Формат вызова: ./ion4 taskFileName outputFileName\n";
            return -1;
        }
        else
        {
            taskFileName = std::string(argv[1]);
            outputFileName = std::string(argv[2]);
        }

        taskType task;
        SimpleParser parser(taskFileName);
        parser.GetTask(task);

        std::vector<unsigned int> Z;
        std::vector<double> x;
        double lgRhoMin, lgRhoMax, lgRhoStep;
        double lgTMin, lgTMax, lgTStep;
        double lgtDiffTmin, lgtDiffTmax;
        double rCoeff;

        if(task["Z"].empty())
        {
            std::cout << "Задайте компоненты смеси\n";
            return -1;
        }
        else
        {
            for(double item : task["Z"]) Z.push_back(item + 0.5);
        }

        if(task["Z"].size() != task["x"].size())
        {
            std::cout << "Задайте концентрации для каждой компоненты\n";
            return -1;
        }
        else
        {
            x = task["x"];
        }

        if(task["lgRho"].size() != 3)
        {
            std::cout << "Задайте диапазон по плотности\n";
            return -1;
        }
        else
        {
            lgRhoMin = task["lgRho"][0];
            lgRhoMax = task["lgRho"][1];
            lgRhoStep = task["lgRho"][2];
        }

        if(task["lgT"].size() != 3)
        {
            std::cout << "Задайте диапазон по ионной температуре\n";
            return -1;
        }
        else
        {
            lgTMin = task["lgT"][0];
            lgTMax = task["lgT"][1];
            lgTStep = task["lgT"][2];
        }

        if(task["lgtDiffT"].size() != 2)
        {
            std::cout << "Задайте диапазон отличия электронной температуры от ионной\n";
            return -1;
        }
        else
        {
            lgtDiffTmin = task["lgtDiffT"][0];
            lgtDiffTmax = task["lgtDiffT"][1];
        }

        if(task["rCoeff"].size() != 1)
        {
            std::cout << "Задайте коэффициент объема ионных остовов\n";
            return -1;
        }
        else
        {
            rCoeff = task["rCoeff"][0];
        }

        if(task["GPaKJ"].size() > 0)
        {
            MixData::SetOutputFormat(MixData::GPaKJ);
        }

        if(task["TPaMJ"].size() > 0)
        {
            MixData::SetOutputFormat(MixData::TPaMJ);
        }

        for(auto &taskItem : task)
        {
            size_t uPos = taskItem.first.find_first_of('_');
            if(uPos != std::string::npos)
            {
                std::string key = taskItem.first.substr(0, uPos);
                std::stringstream ss(taskItem.first.substr(uPos + 1));
                unsigned int z;ss >> z;
                if(key == "extFi")
                {
                    saha::SetExternalFi(z, taskItem.second);
                }
                else if(key == "extG")
                {
                    saha::SetExternalG(z, taskItem.second);
                }
            }
        }

        calculatorMix(Z, x, rCoeff, lgRhoMin, lgRhoMax, lgRhoStep, lgTMin, lgTMax, lgTStep, lgtDiffTmin, lgtDiffTmax, outputFileName);
	}
	catch (std::exception& r)
	{
		printf("\n%s", r.what());
	}
	
    return 0;
}
