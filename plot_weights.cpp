/*
c++ -o plot_weights plot_weights.cpp `root-config --glibs --cflags`
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
		
	TApplication* myapp = new TApplication ("myapp", NULL, NULL);
	TCanvas* cnv = new TCanvas("cnv","cnv",0,0,1200,400);
	cnv->Divide(3,1);

	string name_files[] = {"ntuple_SMlimit.root","ntuple_RcW_0p3.root","ntuple_RcW_0p3.root"};
	string name_ntuples[] = {"SSeu_SMlimit","SSeu_RcW_bsm_0p3","SSeu_RcW_int_0p3"};
	string name_global_numbers[] = {"SSeu_SMlimit_nums","SSeu_RcW_bsm_0p3_nums","SSeu_RcW_int_0p3_nums"};
	string name_histograms[] = {"SM simulation","BSM (quadratic term)","BSM (interference term)"};
	
	for (int j = 0; j < 3; j++) // j = 0,1,2: plots of SM simulation, BSM (quadratic term), BSM (interference term)
	{
		gStyle->SetHistFillColor(j+3);
		TFile* myfile = new TFile(name_files[j].c_str());
		TTreeReader reader (name_ntuples[j].c_str(), myfile);
		TTreeReaderValue<float> var1 (reader, kinetic_variable);
		TTreeReaderValue<float> var2 (reader, "w"); //weights branch
		vector<float> values;
		vector<float> weights;
		
		while (reader.Next ()) 
		{	
			values.push_back(*var1);
			weights.push_back(*var2);
		}		
		
		double min = *min_element(values.begin(), values.end());
		double max = *max_element(values.begin(), values.end());
		TH1F* histo = new TH1F (kinetic_variable, name_histograms[j].c_str(), 100, min, max);
		TH1F* global_numbers = (TH1F*) myfile->Get(name_global_numbers[j].c_str()) ;
		double cross_section = global_numbers->GetBinContent(1);
		double sum_weights_total = global_numbers->GetBinContent(2);
		double sum_weights_selected = global_numbers->GetBinContent(3);
		double luminosity = 100;
		double normalization = cross_section*luminosity/sum_weights_total;	

		for (int i = 0; i < values.size(); i++) 
		{
			histo->Fill(values[i],weights[i]);
		}

		histo->Scale(normalization);
		
		cnv->cd(j+1);
		histo->Draw("HIST");
		histo->GetXaxis()->SetTitle(kinetic_variable);
		cnv->Modified();
		cnv->Update();
	
	}

	myapp->Run();

	return 0;

}

