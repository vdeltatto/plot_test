/*
c++ -o plot_tot plot_tot.cpp `root-config --glibs --cflags`
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
	TCanvas* cnv = new TCanvas("cnv","cnv",0,0,1200,400);
	cnv->Divide(3,1);

	string wilson_values[] = {"0p1","0p3","1"};
	string titles[] = {"cW = 0.1", "cW = 0.3", "cW = 1"};
	string words[] = {"ntuple_RcW_",".root","SSeu_RcW_bsm_","SSeu_RcW_int_","_nums"};
	
	vector<TH1F*> histos;
	float min_tot, max_tot, bin_width;
	
	for (int k = 0; k < 3; k++) // k = 0,1,2: cW = 0.1, 0.3, 1.
	{
		string name_files[] = {"ntuple_SMlimit.root", words[0] + wilson_values[k] + words[1],
			words[0] + wilson_values[k] + words[1]};
		string name_ntuples[] = {"SSeu_SMlimit", words[2] + wilson_values[k], 
			words[3] + wilson_values[k]};
		string name_global_numbers[] = {"SSeu_SMlimit_nums", words[2] + wilson_values[k] + words[4],
			words[3] + wilson_values[k] + words[4]};

		vector<float> values[3]; // to contain data of SM simulation, BSM (quadratic term), BSM (interference term)
		vector<float> weights[3];
		
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
			
			float min = *min_element(values[j].begin(), values[j].end());
			float max = *max_element(values[j].begin(), values[j].end());
			if (j == 0)
			{
				min_tot = min; 
				max_tot = max;
			}
			if (min < min_tot) min_tot = min;
			if (max > max_tot) max_tot = max;	

			myfile->Close();
		}
		for (int j = 0; j < 3; j++) // j = 0,1,2: SM simulation, BSM (quadratic term), BSM (interference term)
		{
			TFile* myfile = new TFile(name_files[j].c_str());
			int Nbins = 100;
			
			if (k == 0 && j == 0) 
			{
				bin_width = (max_tot - min_tot)/Nbins;
			}
			else if (k == 1 || k == 2)
			{
				Nbins = floor((max_tot - min_tot)/bin_width);
			}

			TH1F* histo = new TH1F ("histo", "histo", Nbins, min_tot, max_tot);
				
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
		histo_sum->SetTitle(titles[k].c_str());
		histo_sum->SetName(kinetic_variable);
		
		cnv->cd(k+1);
		
		histo_sum->Draw("HIST");
		histo_sum->GetXaxis()->SetTitle(kinetic_variable);
		histo_sum->GetYaxis()->SetTitle("# events");
		histos[0]->Draw("HIST SAME");
		histo_sum->SetLineColor(kRed);
		histo_sum->SetLineWidth(2.);
		TLegend *legend = new TLegend(0.2,0.8,0.65,0.9); 			
		legend->AddEntry("histo","SM","l");
		legend->AddEntry(kinetic_variable,"SM + BSM + interference","l");
		legend->SetTextFont(42);
	    	legend->SetTextSize(0.03);
		legend->Draw();

		cnv->Modified();
		cnv->Update();	

		histos.clear();
	}

	myapp->Run();

	return 0;

}

