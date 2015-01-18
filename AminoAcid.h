#ifndef AminoAcid_H
#define AminoAcid_H
#include <string>

using namespace std;

class AminoAcid{
    
    private: 
        string name, abbrev, formula;
        
    public: 
        AminoAcid(string abbrev, string name, string formula);
        string getName();
        void setName(string thisName);
        string getAbbrev();
        void setAbbrev(string thisAbbrev);
        string getFormula();
        void setFormula(string thisFormula);
};

#endif
