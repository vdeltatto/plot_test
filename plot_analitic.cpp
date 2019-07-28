/*
c++ -o plot_analitic plot_analitic.cpp `root-config --glibs --cflags`
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

	string titles[] = {"cW = 0.1", "cW = 0.3 (with scaling relations from cW = 0.1)", "cW = 1 (with scaling relations from cW = 0.1)"};
	string name_files[] = {"ntuple_SMlimit.root", "ntuple_RcW_0p1.root","ntuple_RcW_0p1.root"};
	string name_ntuples[] = {"SSeu_SMlimit","SSeu_RcW_bsm_0p1","SSeu_RcW_int_0p1"};
	string name_global_numbers[] = {"SSeu_SMlimit_nums", "SSeu_RcW_bsm_0p1_nums","SSeu_RcW_int_0p1_nums"};

	vector<TH1F*> histos;
	float min_tot, max_tot;

	vector<float> values[3];
	vector<float> weights[3];
		
	for (int j = 0; j < 3; j++) // SM simulation, BSM (quadratic term), BSM (interference term) for cW = 0.1
	{
		TFile* myfile = new TFile(name_files[j].c_str());
		TTreeReader reader (name_ntuples[j].c_str(), myfile);
		TTreeReaderValue<float> var1 (reader, kinetic_variable);
		TTreeReaderValue<float> var2 (reader, "w"); //weights branch
		
		while (reader.Next ()) 
		{	
			values[j].push_back(*var1);
			weights[j].push_back(*var2);
		}		
		
		double min = *min_element(values[j].begin(), values[j].end());
		float max = *max_element(values[j].begin(), values[j].end());
		if (j == 0)
		{
			min_tot = min; 
			max_tot = max;
		}
		else if (min < min_tot) min_tot = min;
		else if (max > max_tot) max_tot = max;	
	}
	for (int j = 0; j < 3; j++) // SM simulation, BSM (quadratic term), BSM (interference term) for cW = 0.1
	{
		TFile* myfile = new TFile(name_files[j].c_str());
		TH1F* histo = new TH1F ("histo", "histo", 100, min_tot, max_tot);
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

	for (int j = 0; j < 3; j++) // j = 0 -> cW = 0.1, j = 1 -> cW = 0.3, j = 2 -> cW = 1
	{ 
		if (j == 1) 
		{
			histos[1]->Scale(0.3*0.3/(0.1*0.1)); // quadratic scaling relation
			histos[2]->Scale(0.3/0.1);           // linear scaling relation
		}
		else if (j == 2)
		{
			histos[1]->Scale(1*1/(0.3*0.3)); // quadratic scaling relation
			histos[2]->Scale(1/0.3);	 // linear scaling relation
		}
	
		TH1F* histo_sum = new TH1F(*histos[0]);
		histo_sum->Add(histos[1]); 
		histo_sum->Add(histos[2]);
		histo_sum->SetTitle(titles[j].c_str());
		histo_sum->SetName(kinetic_variable);
		
		cnv->cd(j+1);
		
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
	}

	myapp->Run();

	return 0;

}

