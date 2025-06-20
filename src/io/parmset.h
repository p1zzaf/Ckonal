#ifndef PARMSET_H
#define PARMSET_H

#include <iostream>
#include <string>
#include <map>
#include <vector>

using namespace std;

class ParmSet
{
    friend std::ostream& operator<< (std::ostream &, const ParmSet&);
    friend std::istream& operator>> (std::istream &, ParmSet&);
    public:
        ParmSet();
        ParmSet(std::string);
        virtual ~ParmSet();
        bool isDef(std::string);

        std::string getString(std::string kw);
        int getInt(std::string kw);
        double getFloat(std::string kw);
        vector<double> getFloats(std::string kw);

        bool getString(std::string kw, std::string &out);
        bool getInt(std::string kw, int &out);
        bool getFloat(std::string kw, double &out);

	void putKeyValue(std::string kw, std::string vl);

	void changeKeyValue(std::string kw, std::string vl);

    protected:
    private:
        std::map<std::string,std::string> table;
};



#endif // PARMSET_H
