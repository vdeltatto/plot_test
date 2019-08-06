/*
Plot of 'kinetic_variable' distribution for cW = 0.05, 0.1, 0.3, 0.4, 1. The distribution 
with cW = 0.3 is plotted with data stored in ntuples, the others are created using scaling relations 
(quadratic for the pure BSM term and linear for the interference one).

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
		
	TApplication* myapp = new TApplication ("myapp", NULL, NULL);
	TCanvas* cnv = new TCanvas("cnv","cnv",0,0,1000,400);
	TCanvas* cnv_logy = new TCanvas("cnv_logy","cnv_logy",0,450,1000,400);
	TCanvas* cnv2 = new TCanvas("cnv2","cnv2",0,0,1000,400);
	TCanvas* cnv2_logy = new TCanvas("cnv2_logy","cnv2_logy",0,450,1000,400);
		
	cnv->Divide(2,1);
	cnv_logy->Divide(2,1);
	cnv2->Divide(2,1);
	cnv2_logy->Divide(2,1);

	string titles[] = {"cW = 0.3", "cW = 0.05 (scaling relations)", "cW = 0.1 (scaling relations)", 
		"cW = 0.4 (scaling relations)", "cW = 1 (scaling relations)"};
	string name_files[] = {"ntuple_SMlimit_HS.root", "ntuple_RcW_0p3_HS.root", "ntuple_RcW_0p3_HS.root"};
	string name_ntuples[] = {"SSeu_SMlimit","SSeu_RcW_bsm_0p3","SSeu_RcW_int_0p3"};
	string name_global_numbers[] = {"SSeu_SMlimit_nums", "SSeu_RcW_bsm_0p3_nums","SSeu_RcW_int_0p3_nums"};
	string name_histograms[] = {"SM", "BSM", "interference"};
	string kinetic_variables[] = {"met","mjj","mll","ptl1","ptl2"};

	vector<TH1F*> histos;
	vector<float> values[3];
	vector<float> weights[3];

	float max_tot;
	float maxima[] = {800, 6000, 1500, 1000, 400};
	for (int i = 0; i < 5; i++) 	
	{
		if (kinetic_variable == kinetic_variables[i])
		{
			max_tot = maxima[i];	
			break;
		}
	}	
		
	for (int j = 0; j < 3; j++) // SM simulation, BSM (quadratic term), BSM (interference term)
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

		int Nbins = 70;
		TH1F* histo = new TH1F ("histo", name_histograms[j].c_str(), Nbins, 0., max_tot);
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

	for (int k = 1; k < 5; k++) // cW = 0.05, 0.1, 0.3, 0.4, 1
	{ 
		THStack* h_stack = new THStack("hs","");
		
		if (k == 1) 
		{
			histos[1]->Scale(0.05*0.05/(0.3*0.3)); // quadratic scaling relation
			histos[2]->Scale(0.05/0.3);           // linear scaling relation
		}
		else if (k == 2)
		{
			histos[1]->Scale(0.1*0.1/(0.05*0.05)); // quadratic scaling relation
			histos[2]->Scale(0.1/0.05);	 // linear scaling relation
		}
		else if (k == 3)
		{
			histos[1]->Scale(0.4*0.4/(0.1*0.1)); // quadratic scaling relation
			histos[2]->Scale(0.4/0.1);	 // linear scaling relation
		}
		else if (k == 4)
		{
			histos[1]->Scale(1/(0.4*0.4)); // quadratic scaling relation
			histos[2]->Scale(1/0.4);	 // linear scaling relation
		}
	
		histos[0]->SetLineColor(kBlue);
		histos[1]->SetLineColor(kRed);
		histos[2]->SetLineColor(kGreen + 1);

		for (int i = 0; i < 3; i++)
		{
			h_stack->Add(histos[i]);
		}

		TH1F* histo_sum = new TH1F(*histos[0]);
		histo_sum->Add(histos[1]); 
		histo_sum->Add(histos[2]);
		histo_sum->SetTitle("SM + BSM + interference");
		histo_sum->SetLineColor(kBlack);
		//histo_sum->SetLineWidth(2);
		h_stack->Add(histo_sum);
		
		if (k < 3) cnv->cd(k);
		else cnv2->cd(k-2);

		h_stack->Draw("HIST NOSTACK");
		TText* T = new TText(); 
		T->SetTextFont(42); 
		T->SetTextAlign(21);
		T->DrawTextNDC(.5,.95,titles[k].c_str());
		h_stack->GetXaxis()->SetTitle(kinetic_variable);
		h_stack->GetYaxis()->SetTitle("# events"); 
		gPad->BuildLegend(0.40,0.70,0.90,0.90,"");
			
		if (k < 3) 
		{
			cnv->Modified();
			cnv->Update();
		}
		else	
		{
			cnv2->Modified();
			cnv2->Update();
		}

		if (k < 3) cnv_logy->cd(k); //logarithmic plot
		else cnv2_logy->cd(k-2);
		
		h_stack->Draw("HIST NOSTACK");
		TText* T_logy = new TText(); 
		T_logy->SetTextFont(42); 
		T_logy->SetTextAlign(21);
		string title = titles[k] + " (logarithmic scale)";
		T_logy->DrawTextNDC(.5,.95,title.c_str());
		h_stack->GetXaxis()->SetTitle(kinetic_variable);
		h_stack->GetYaxis()->SetTitle("# events"); 
		gPad->BuildLegend(0.40,0.70,0.90,0.90,"");
		gPad->SetLogy();

		if (k < 3) 
		{
			cnv_logy->Modified();
			cnv_logy->Update();
		}
		else	
		{
			cnv2_logy->Modified();
			cnv2_logy->Update();
		}
	
		histos.clear();
	}
	//To save the plots
	/*string name1_png = string(kinetic_variable) + "_1_analitic.png";
	string name1_logy_png = string(kinetic_variable) + "_1_analitic_log.png";
	string name2_png = string(kinetic_variable) + "_2_analitic.png";
	string name2_logy_png = string(kinetic_variable) + "_2_analitic_log.png";
	cnv->Print(name1_png.c_str(), "png");
	cnv_logy->Print(name1_logy_png.c_str(), "png");
	cnv2->Print(name2_png.c_str(), "png");
	cnv2_logy->Print(name2_logy_png.c_str(), "png");*/

	myapp->Run();

	return 0;

}

