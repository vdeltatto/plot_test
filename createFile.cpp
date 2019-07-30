/*
c++ -o createFile createFile.cpp `root-config --glibs --cflags` 
*/

#include <iostream>
#include <string>

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


int main () 
{
	TFile f ("histo_0p3.root","recreate");
	const char* kinetic_variable; //possible variables: met, mjj, mll, ptl1, ptl2

	string name_histograms[] = {"sm_", "quadratic_", "linear_"};
	string kinetic_variables[] = {"met","mjj","mll","ptl1","ptl2"};
	string name_files[] = {"ntuple_SMlimit.root", "ntuple_RcW_0p3.root", "ntuple_RcW_0p3.root"};
	string name_ntuples[] = {"SSeu_SMlimit", "SSeu_RcW_bsm_0p3", "SSeu_RcW_int_0p3"};
	string name_global_numbers[] = {"SSeu_SMlimit_nums", "SSeu_RcW_bsm_0p3_nums", "SSeu_RcW_int_0p3_nums"};

	float min, max;
	vector<float> values; 
	vector<float> weights;
	int Nbins = 100;
	vector<TH1F*> histos;
		
	for (int k = 0; k < 5; k++) 
	{
		kinetic_variable = kinetic_variables[k].c_str();
		for (int j = 0; j < 3; j++) // j = 0,1,2: SM simulation, BSM (quadratic term), BSM (interference term)
		{ 
			TFile* myfile = new TFile(name_files[j].c_str());
			TTreeReader reader (name_ntuples[j].c_str(), myfile);
			TTreeReaderValue<float> var1 (reader, kinetic_variable);
			TTreeReaderValue<float> var2 (reader, "w"); //weights branch
			
			while (reader.Next()) 
			{	
				values.push_back(*var1);
				weights.push_back(*var2);
			}		
			
			min = *min_element(values.begin(), values.end());
			max = *max_element(values.begin(), values.end());

			string name = name_histograms[j] + string(kinetic_variable);
			TH1F* histo = new TH1F (name.c_str(), kinetic_variable, Nbins, min, max);	
				
			TH1F* global_numbers = (TH1F*) myfile->Get(name_global_numbers[j].c_str()) ;
			float cross_section = global_numbers->GetBinContent(1);
			float sum_weights_total = global_numbers->GetBinContent(2);
			float sum_weights_selected = global_numbers->GetBinContent(3);
			float luminosity = 100;
			float normalization = cross_section*luminosity/sum_weights_total;

			f.cd();

			for (int i = 0; i < values.size(); i++) 
			{
				histo->Fill(values[i],weights[i]);
			}

			histo->Scale(normalization);
			histo->Write(name.c_str());
			f.Write();

			values.clear();
			weights.clear();	
		}
	}

	return 0 ;
}

