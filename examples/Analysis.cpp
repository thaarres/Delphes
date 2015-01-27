#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"

#include "TH2.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

#include <iostream>
#include <utility>
#include <vector>
#include <map>


#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

using namespace std;

//------------------------------------------------------------------------------

// Analysis macro
#include "Analysis.C"

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char *appName = "Analysis";

  if(argc != 3)
  {
    cout << " Usage : " << appName << " input_file output_file" << endl;
    cout << " Usage : " << appName << " input_file output_file" << endl;
    cout << " input_file - input file in ROOT format ('Delphes' tree)" << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  TString inputFile(argv[1]);
  TString outputFile(argv[2]);

//------------------------------------------------------------------------------

// Call the macros main function (in Analysis.C)

  Analysis(inputFile, outputFile);

//------------------------------------------------------------------------------

}


