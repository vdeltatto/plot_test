/*

Returns the RMS of the total distribution (SM + BSM + interference)

c++ -o getRMS getRMS.cpp `root-config --glibs --cflags`

./getRMS mll

*/

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <TFile.h>
#include <TNtuple.h>
#include <TTreeReader.h>
#include <TH1.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <THStack.h>
#include <TText.h>

using namespace std ;


int main (int argc, char** argv) 
{
	const char* kinetic_variable; //possible variables: met, mjj, mll, ptl1, ptl2
	if (argc == 1) 
	{
		kinetic_variable = "met";
	}
	else 
	{
		kinetic_variable = argv[1];
	}
		
	TH1::SetDefaultSumw2();

	string wilson_values[] = {"0p05","0p1","0p3","0p4","1"};
	string titles[] = {"cW = 0.05", "cW = 0.1", "cW = 0.3 (high statistics)","cW = 0.4", "cW = 1"};
	string words[] = {"ntuple_RcW_",".root","SSeu_RcW_bsm_","SSeu_RcW_int_","_nums","_HS.root"};
	string name_histograms[] = {"SM", "BSM", "interference"};
	string kinetic_variables[] = {"met","mjj","mll","ptl1","ptl2"};
	
	vector<TH1F*> histos;
	
	float max_tot[5];
	
	for (int k = 0; k < 5; k++) // cW = 0.05, 0.1, 0.3, 0.4, 1.
	{
	string name_files[3];
		if (k != 2) //for cW = 0.3 high statistics (HS)
		{
			name_files[0] = "ntuple_SMlimit_HS.root";
			name_files[1] = words[0] + wilson_values[k] + words[1];
			name_files[2] = words[0] + wilson_values[k] + words[1];
		}
		else 
		{
			name_files[0] = "ntuple_SMlimit_HS.root";
			name_files[1] = words[0] + wilson_values[k] + words[5];
			name_files[2] = words[0] + wilson_values[k] + words[5];
		}
		string name_ntuples[] = {"SSeu_SMlimit", words[2] + wilson_values[k], 
			words[3] + wilson_values[k]};
		string name_global_numbers[] = {"SSeu_SMlimit_nums", words[2] + wilson_values[k] + words[4],
			words[3] + wilson_values[k] + words[4]};

		vector<float> values[3]; // to contain data of SM simulation, BSM (quadratic term), BSM (interference term)
		vector<float> weights[3];
		int Nbins = 200;
		float max = 0.;
		
		for (int j = 0; j < 3; j++) // j = 0,1,2: SM simulation, BSM (quadratic term), BSM (interference term)
		{
			TFile* myfile = new TFile(name_files[j].c_str());
			TTreeReader reader (name_ntuples[j].c_str(), myfile);
			TTreeReaderValue<float> var1 (reader, kinetic_variable);
			TTreeReaderValue<float> var2 (reader, "w"); //weights branch
			
			while (reader.Next()) 
			{	
				values[j].push_back(*var1);
				weights[j].push_back(*var2);
			}		

			if (max < *max_element(values[j].begin(), values[j].end())) max = *max_element(values[j].begin(), values[j].end());

		}
		for (int j = 0; j < 3; j++) // j = 0,1,2: SM simulation, BSM (quadratic term), BSM (interference term)
		{
			TFile* myfile = new TFile(name_files[j].c_str());
			TTreeReader reader (name_ntuples[j].c_str(), myfile);
			
			TH1F* histo = new TH1F ("histo", name_histograms[j].c_str(), Nbins, 0., max);		
				
			TH1F* global_numbers = (TH1F*) myfile->Get(name_global_numbers[j].c_str()) ;
			float cross_section = global_numbers->GetBinContent(1);
			float sum_weights_total = global_numbers->GetBinContent(2);
			float sum_weights_selected = global_numbers->GetBinContent(3);
			float luminosity = 100;
			float normalization = cross_section*luminosity/sum_weights_total;	

			for (int i = 0; i < values[j].size(); i++) 
			{
				histo->Fill(values[j][i],weights[j][i]);
			}

			histo->Scale(normalization);

			histos.push_back(histo);	

			values[j].clear();
			weights[j].clear();	
		}
		
		TH1F* histo_sum = new TH1F(*histos[0]);
		histo_sum->Add(histos[1]);
		histo_sum->Add(histos[2]);

		float RMS = histo_sum->GetRMS();
		cout << "RMS (" << titles[k] << "): " << RMS << endl;
	}

	return 0;

}
