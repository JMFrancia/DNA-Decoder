#include "AminoAcid.h"

using namespace std;

AminoAcid::AminoAcid(string abbrev, string name, string formula){
    setName(name);
    setAbbrev(abbrev);
    setFormula(formula);    
}

void AminoAcid::setName(string thisName){
    name = thisName;    
}

void AminoAcid::setAbbrev(string thisAbbrev){
    abbrev = thisAbbrev;    
}

void AminoAcid::setFormula(string thisFormula){
    formula = thisFormula;
}

string AminoAcid::getName(){return name;}
string AminoAcid::getAbbrev(){return abbrev;}
string AminoAcid::getFormula(){return formula;}
