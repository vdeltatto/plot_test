/*
c++ -o compare compare.cpp `root-config --glibs --cflags`
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
	TCanvas* cnv = new TCanvas("cnv","cnv",0,0,1000,400);
	cnv->Divide(2,1);

	string wilson_values[] = {"0p1","0p3","1"};
	string titles[] = {"cW = 0.1", "cW = 0.3 (data - expected_values_from_scaling_relations)", 
		"cW = 1 (data - expected_values_from_scaling_relations)"};
	string words[] = {"ntuple_RcW_",".root","SSeu_RcW_bsm_","SSeu_RcW_int_","_nums"};
	string kinetic_variables[] = {"met","mjj","mll","ptl1","ptl2"};
	
	vector<TH1F*> histos, histos_analitic;
	float max;
	float maxima[] = {3500, 8500, 7000, 4000, 3500};
	for (int i = 0; i < 5; i++) 	
	{
		if (kinetic_variable == kinetic_variables[i])
		{
			max = maxima[i];	
			break;
		}
	}	

	for (int k = 0; k < 3; k++) // cW = 0.1, cW = 0.3, cW = 1
	{
		string name_files[] = {"ntuple_SMlimit.root", words[0] + wilson_values[k] + words[1],
			words[0] + wilson_values[k] + words[1]};
		string name_ntuples[] = {"SSeu_SMlimit", words[2] + wilson_values[k], 
			words[3] + wilson_values[k]};
		string name_global_numbers[] = {"SSeu_SMlimit_nums", words[2] + wilson_values[k] + words[4],
			words[3] + wilson_values[k] + words[4]};
		vector<float> values;
			vector<float> weights;
		
		
		for (int j = 0; j < 3; j++) // j = 0,1,2: SM simulation, BSM (quadratic term), BSM (interference term)
		{
			
			TFile* myfile = new TFile(name_files[j].c_str());
			TTreeReader reader (name_ntuples[j].c_str(), myfile);
			TTreeReaderValue<float> var1 (reader, kinetic_variable);
			TTreeReaderValue<float> var2 (reader, "w"); //weights branch
			
			while (reader.Next ()) 
			{	
				values.push_back(*var1);
				weights.push_back(*var2);
			}		
			
			int Nbins = 100;

			TH1F* histo = new TH1F ("histo", "histo", Nbins, 0, max);
			
			TH1F* global_numbers = (TH1F*) myfile->Get(name_global_numbers[j].c_str());
			float cross_section = global_numbers->GetBinContent(1);
			float sum_weights_total = global_numbers->GetBinContent(2);
			float sum_weights_selected = global_numbers->GetBinContent(3);
			float luminosity = 100;
			float normalization = cross_section*luminosity/sum_weights_total;	

			for (int i = 0; i < values.size(); i++) 
			{
				histo->Fill(values[i],weights[i]);
			}

			histo->Scale(normalization);
			
			if (k == 0) histos_analitic.push_back(histo); 	
			else histos.push_back(histo);					
					
			values.clear();		
			weights.clear();
		}

		if (k == 1 || k == 2) 
		{
			if (k == 1)
			{
				histos_analitic[1]->Scale(0.3*0.3/(0.1*0.1)); // quadratic scaling relation (pure BSM term)
				histos_analitic[2]->Scale(0.3/0.1);           // linear scaling relation (interference term)
			}
			else if (k == 2)
			{
				histos_analitic[1]->Scale(1./(0.3*0.3)); // quadratic scaling relation (pure BSM term)
				histos_analitic[2]->Scale(1/0.3);        // linear scaling relation (interference term)   
			}

			TH1F* histo_sum = new TH1F(*histos[0]);
			histo_sum->Add(histos[1]); 
			histo_sum->Add(histos[2]);
			histo_sum->SetTitle(titles[k].c_str());
			histo_sum->SetName(kinetic_variable);

			TH1F* histo_sum_analitic = new TH1F(*histos_analitic[0]);
			histo_sum_analitic->Add(histos_analitic[1]); 
			histo_sum_analitic->Add(histos_analitic[2]);
			histo_sum_analitic->SetTitle(titles[k].c_str());
			histo_sum_analitic->SetName(kinetic_variable);

			TH1F* compare = new TH1F(*histo_sum - *histo_sum_analitic);
			
			cnv->cd(k);
		
			compare->Draw();
			compare->GetXaxis()->SetTitle(kinetic_variable);
			compare->GetYaxis()->SetTitle("difference between events");

			cnv->Modified();
			cnv->Update();
		}
		histos.clear();
	}

	string name_png = string(kinetic_variable) + "_compare.png";
	cnv->Print(name_png.c_str(),"png");

	myapp->Run();

	return 0;

}

