/*
J.M.Francia
CS135
HW Assignment #2: Amino Acids

Run HW_Amino with a file name as a parameter: ./HW_Amino "thisFile.txt". 
Or, alternatively, use "test" as a parameter to use the DNA strand from the assignment sheet

This program takes a text file as an input, containing a string of character to make a DNA strand.
It will transcribe the DNA to RNA, then translate the RNA to display the sequence of amino acids used.
It will then provide the name and formula for each amino acid used.

*/


#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "AminoAcid.h"

using namespace std;

string getDNAFile(string filename); //Retrieves DNA strand from text file
string getTestDNA(); //retrieves test DNA strand from assignment sheet
string transcribeRNA(string DNA); //Transcribes DNA strand, returns as RNA strand
string codonToAA(string codon, string codonDB[][2]); //Accepts codon, returns abbreviation of AA
AminoAcid * getAA(string abbrev, AminoAcid * aaPointers[]); //Accepts the abbreviation of an AA, returns a pointer to that AA
bool isStartCodon(string codon); //Will check to see if a codon is a start codon
bool isStopCodon(string codon);  //Will check to see if a codon is a stop codon
bool isOnList(string AA, string totalAA[], int totalCounter); //Will check to see if string AA exists on list totalAA

const int TOTAL_CODONS = 64;

int main(int argc, char *argv[])
{    
    
    //Test to check if number of arguments is correct
    if(argc != 2){
        cout << "\nTo use this program, type ./HW_Amino <filename>" << endl;
        cout <<"Alternatively, type ./HW_Amino \"test\" to run program with sample DNA strand from assignment sheet" << endl;
        return 1;
    }
    
    
    //Codon database
    string codonDB[TOTAL_CODONS][2] = {
                             {"uuu","Phe"},{"uuc","Phe"},
                             {"uua","Leu"},{"uug","Leu"},{"cuu","Leu"},{"cuc","Leu"},{"cua","Leu"},{"cug","Leu"},
                             {"auu","Ile"},{"auc","Ile"},{"aua","Ile"},
                             {"aug","Met"},                           
                             {"guu","Val"},{"guc","Val"},{"gua","Val"},{"gug","Val"},
                             {"ucu","Ser"},{"ucc","Ser"},{"uca","Ser"},{"ucg","Ser"},{"agu","Ser"},{"agc","Ser"},
                             {"ccu","Pro"},{"ccc","Pro"},{"cca","Pro"},{"ccg","Pro"},
                             {"acu","Thr"},{"acc","Thr"},{"aca","Thr"},{"acg","Thr"},
                             {"gcu","Ala"},{"gcc","Ala"},{"gca","Ala"},{"gcg","Ala"},
                             {"uau","Tyr"},{"uac","Tyr"},
                             {"uaa","TER"},{"uag","TER"},{"uga","TER"},
                             {"cau","His"},{"cac","His"},
                             {"caa","Gln"},{"cag","Gln"},
                             {"aau","Asn"},{"aac","Asn"},
                             {"aaa","Lys"},{"aag","Lys"},
                             {"gau","Asp"},{"gac","Asp"},
                             {"gaa","Glu"},{"gag","Glu"},
                             {"ugu","Cys"},{"ugc","Cys"},
                             {"ugg","Trp"},
                             {"cgu","Arg"},{"cgc","Arg"},{"cga","Arg"},{"cgg","Arg"},{"aga","Arg"},{"agg","Arg"},
                             {"ggu","Gly"},{"ggc","Gly"},{"gga","Gly"},{"ggg","Gly"}
                            };
    
    //Amino Acids   
    AminoAcid Phe("Phe","Phenylalanine","C6H5CH2CH(NH2)COOH"),
              Leu("Leu","Leucine","HO2CCH(NH2)CH2CH(CH3)2"),
              Ile("Ile","Isoleucine","HO2CCH(NH2)CH(CH3)CH2CH3"),
              Met("Met", "Methionine", "HO2CCH(NH2)CH2CH2SCH3"),
              Val("Val", "Valine", "C5H11NO2"),
              Ser("Ser", "Serine", "HO2CCH(NH2)CH2OH"),
              Pro("Pro", "Proline", "C5H9NO2"),
              Thr("Thr", "Threonine", "HO2CCH(NH2)CH(OH)CH3"),
              Ala("Ala", "Alanine", "CH3CH(NH2)COOH"),
              Tyr("Tyr", "Tyrosine", "C9H11NO3"),
              His("His", "Histidine", "C6H9N3O2"),
              Gln("Gln", "Glycine", "C2H5NO2"),
              Asn("Asn", "Asparagine", "C4H8N2O3"),
              Lys("Lys", "Lysine", "HO2CCH(NH2)(CH2)4NH2"),
              Asp("Asp", "Aspartic Acid", "HOOCCH(NH2)CH2COOH"),
              Glu("Glu", "Glutamic Acid", "C5H9NO4"),
              Cys("Cys", "Cysteine", "HO2CCH(NH2)CH2SH"),
              Trp("Trp", "Tryptophan", "C11H12N2O2"),
              Arg("Arg", "Arginine", "C6H14N4O2"),
              Gly("Gly", "Glycine", "NH2CH2COOH");
    
    //Amino Acid pointers
    AminoAcid * Phe_P = &Phe, * Leu_P = &Leu, * Ile_P = &Ile, * Met_P = &Met,
              * Val_P = &Val, * Ser_P = &Ser, * Pro_P = &Pro, * Thr_P = &Thr,
              * Ala_P = &Ala, * Tyr_P = &Tyr, * His_P = &His, * Gln_P = &Gln,
              * Asn_P = &Asn, * Lys_P = &Lys, * Asp_P = &Asp, * Glu_P = &Glu,
              * Cys_P = &Cys, * Trp_P = &Trp, * Arg_P = &Arg, * Gly_P = &Gly;
    
    //Array of all AA pointers so functions can use it
    AminoAcid * aaPointers[] = {Phe_P, Leu_P, Ile_P, Met_P, Val_P, Ser_P, Pro_P,
                                Thr_P, Ala_P, Tyr_P, His_P, Gln_P, Asn_P, Lys_P,
                                Asp_P, Glu_P, Cys_P, Trp_P, Arg_P, Gly_P};
 
    string DNA, RNA, thisCodon, thisAcid, totalAcids[20], mainInput = argv[1];
    int n = 0, totalAcidCounter = 0; 
    bool recording = false, isFirstInSeq = true;
    AminoAcid * thisAcid_P;

    //Retrieve and transcribe DNA .. if arg is "test", use strand from assignment
    
    if(mainInput.compare("test") == 0) 
        DNA = getTestDNA();
    else
        DNA = getDNAFile(argv[1]);                                    
    
    RNA = transcribeRNA(DNA);
    
    //Print out translated RNA strand
    
    
    cout << endl;
    while(n < RNA.size()){     
        thisCodon = RNA.substr(n,3);
        if(isStartCodon(thisCodon))
            recording = true;
        if(isStopCodon(thisCodon)){
            recording = false;
            isFirstInSeq = true;
            cout << endl;
        } 
        if(recording){
            thisAcid = codonToAA(thisCodon,codonDB);
            if(!isOnList(thisAcid,totalAcids,totalAcidCounter)){
                totalAcids[totalAcidCounter] = thisAcid;
                totalAcidCounter++;
            }
            if(isFirstInSeq){
                cout << thisAcid;
                isFirstInSeq = false;
            } else
                cout << "-" << thisAcid;
            n +=3;
        }
        else 
            n++;
    }
    
    //Print out info for all AA's used in this RNA strand
    
    cout << endl;    
    for(int n = 0; n < totalAcidCounter; n++){
        thisAcid_P = getAA(totalAcids[n],aaPointers);
        cout << (*thisAcid_P).getAbbrev()<< "  " << (*thisAcid_P).getName() << "\t" << (*thisAcid_P).getFormula();
        cout << endl;
    }
                          
    return 0;
}

//FUNCTION DEFINITIONS  ---------------------------------------------------------

bool isOnList(string AA, string totalAA[], int totalCounter){
    bool result = false;
    for(int n = 0; n < totalCounter; n++){
        if(AA.compare(totalAA[n]) == 0)
            result = true;
    }
    return result;
}

string getTestDNA(){
    return "atggtttatggtctctgaattaatctccatgttttatcactaa";
}

string getDNAFile(string filename){
    string dat;
    ifstream inputFile;
    inputFile.open(filename.c_str());
    if(inputFile){
         inputFile >> dat;
         getline(inputFile,dat);
         inputFile.close();
    } else
        cout << "ERROR: Cannot find file " << filename << endl;
    return dat;
}

string transcribeRNA(string DNA){
    for(int n = 0; n < DNA.size(); n++){
        if(DNA[n] == 't')
            DNA[n] = 'u';
    }
    return DNA;
}

bool isStartCodon(string codon){
    bool result = false;
    if(codon.compare("aug") == 0)  
        result = true;
    return result;  
}

bool isStopCodon(string codon){
    bool result = false;
    if(codon.compare("uaa") == 0 || codon.compare("uag") == 0 || codon.compare("uga") == 0)  
        result = true;
    return result;  
}

string codonToAA(string codon, string codonDB[][2]){
    string result;
    bool found = false;
    for(int n = 0; n < TOTAL_CODONS; n++){
        if(codon.compare(codonDB[n][0]) == 0){
            result =  codonDB[n][1];
            found = true;
        }
    }
    if(!found)
        cout << "ERROR: Cannot find AA for codon << " << codon << endl;
    return result;   
}

AminoAcid * getAA(string abbrev, AminoAcid * aaPointers[]){
    AminoAcid * result;
    if(abbrev.compare("Phe") == 0)
        result = aaPointers[0];
    if(abbrev.compare("Leu") == 0)
        result = aaPointers[1];
    if(abbrev.compare("Ile") == 0)
        result = aaPointers[2];
    if(abbrev.compare("Met") == 0)
        result = aaPointers[3];
    if(abbrev.compare("Val") == 0)
        result = aaPointers[4];
    if(abbrev.compare("Ser") == 0)
        result = aaPointers[5];
    if(abbrev.compare("Pro") == 0)
        result = aaPointers[6];
    if(abbrev.compare("Thr") == 0)
        result = aaPointers[7];
    if(abbrev.compare("Ala") == 0)
        result = aaPointers[8];
    if(abbrev.compare("Tyr") == 0)
        result = aaPointers[9];
    if(abbrev.compare("His") == 0)
        result = aaPointers[10];
    if(abbrev.compare("Gln") == 0)
        result = aaPointers[11];
    if(abbrev.compare("Asn") == 0)
        result = aaPointers[12];
    if(abbrev.compare("Lys") == 0)
        result = aaPointers[13];
    if(abbrev.compare("Asp") == 0)
        result = aaPointers[14];
    if(abbrev.compare("Glu") == 0)
        result = aaPointers[15];
    if(abbrev.compare("Cys") == 0)
        result = aaPointers[16];
    if(abbrev.compare("Trp") == 0)
        result = aaPointers[17];
    if(abbrev.compare("Arg") == 0)
        result = aaPointers[18];
    if(abbrev.compare("Gly") == 0)
        result = aaPointers[19];
    return result;
}
