/*
Plot of 'kinetic_variable' distributions for cW = 0.1, 0.3, 1 (with data stored in ntuples)

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
	TCanvas* cnv = new TCanvas("cnv","cnv",0,0,1200,400);
	TCanvas* cnv_logy = new TCanvas("cnv_logy","cnv_logy",0,450,1200,400); 
	TCanvas* zoom = new TCanvas("zoom","zoom",0,0, 1000, 500);
	
	cnv->Divide(3,1);
	cnv_logy->Divide(3,1);

	string wilson_values[] = {"0p1","0p3","1"};
	string titles[] = {"cW = 0.1", "cW = 0.3", "cW = 1"};
	string words[] = {"ntuple_RcW_",".root","SSeu_RcW_bsm_","SSeu_RcW_int_","_nums"};
	string name_histograms[] = {"SM", "BSM", "interference"};
	string kinetic_variables[] = {"met","mjj","mll","ptl1","ptl2"};
	
	vector<TH1F*> histos;
	vector<TH1F*> histos_zoom; //to zoom the critical part of mjj distributon (cW = 1)
	
	float max_tot[3];
	float maxima[][5] = {{1000, 7000, 1500, 1000, 400},{1000, 7000, 1500, 1000, 400},
		{1500, 7000, 3500, 2300, 1200}}; //maxima in the histograms
	for (int i = 0; i < 5; i++) 	
	{
		if (kinetic_variable == kinetic_variables[i])
		{
			for (int n = 0; n < 3; n++) max_tot[n] = maxima[n][i];	
			break;
		}
	}	
	
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
		THStack* h_stack = new THStack("hs","");
		THStack* h_stack_zoom = new THStack("hs_zoom","");
		
		
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
			
			int Nbins = 70;
			
			TH1F* histo = new TH1F ("histo", name_histograms[j].c_str(), Nbins, 0., max_tot[k]);
			TH1F* histo_zoom = new TH1F ("histo_zoom", name_histograms[j].c_str(), 100, 70, 90);			
				
			TH1F* global_numbers = (TH1F*) myfile->Get(name_global_numbers[j].c_str()) ;
			float cross_section = global_numbers->GetBinContent(1);
			float sum_weights_total = global_numbers->GetBinContent(2);
			float sum_weights_selected = global_numbers->GetBinContent(3);
			float luminosity = 100;
			float normalization = cross_section*luminosity/sum_weights_total;	

			for (int i = 0; i < values[j].size(); i++) 
			{
				histo->Fill(values[j][i],weights[j][i]);
				if (k == 2 && kinetic_variable == kinetic_variables[1])
				{
					histo_zoom->Fill(values[j][i],weights[j][i]);
				}
			}

			histo->Scale(normalization);
			histos.push_back(histo);	

			if (k == 2 && kinetic_variable == kinetic_variables[1]) //cW = 1, mjj
			{
				histo_zoom->Scale(normalization);
				histos_zoom.push_back(histo_zoom);	
			}	 

			values[j].clear();
			weights[j].clear();	
		}

		histos[0]->SetLineColor(kBlue);
		histos[1]->SetLineColor(kRed);
		histos[2]->SetLineColor(kGreen +1);
		
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
		
		cnv->cd(k+1);
		
		h_stack->Draw("HIST NOSTACK");
		TText* T = new TText(); 
		T->SetTextFont(42); 
		T->SetTextAlign(21);
		T->DrawTextNDC(.5,.95,titles[k].c_str());
		h_stack->GetXaxis()->SetTitle(kinetic_variable);
		h_stack->GetYaxis()->SetTitle("# events"); 
		gPad->BuildLegend(0.40,0.70,0.90,0.90,"");

		cnv->Modified();
		cnv->Update();	

		cnv_logy->cd(k+1); //logarithmic plot
		
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

		cnv_logy->Modified();
		cnv_logy->Update();	

		histos.clear();		
	
		if (k == 2 && kinetic_variable == kinetic_variables[1]) // zoom in the range with singularity
		{						 	// for cW = 1
			int color_SM  = kGreen - 8 ;
			int color_INT = kOrange - 4 ;
			int color_BSM = kAzure - 9 ;

			histos_zoom[0]->SetFillColor(color_SM);
			histos_zoom[1]->SetFillColor(color_INT);
			histos_zoom[2]->SetFillColor(color_BSM);
			
			for (int i = 0; i < 3; i++)
			{
				h_stack_zoom->Add(histos_zoom[i]);
			}
			
			zoom->Divide(2,1);
			zoom->cd(1);
		
			h_stack_zoom->Draw("HIST NOSTACK");
			TText* T = new TText(); 
			T->SetTextFont(42); 
			T->SetTextAlign(21);
			T->DrawTextNDC(.5,.95,"zoom for cW = 1");
			h_stack_zoom->GetXaxis()->SetTitle(kinetic_variable);
			h_stack_zoom->GetYaxis()->SetTitle("# events"); 
			gPad->BuildLegend(0.10,0.76,0.4,0.90,"");

			zoom->Modified();
			zoom->Update();	

			zoom->cd(2);			
			h_stack_zoom->Draw("HIST NOSTACK");
			TText* T_logy = new TText(); 
			T_logy->SetTextFont(42); 
			T_logy->SetTextAlign(21);
			T_logy->DrawTextNDC(.5,.95,"zoom for cW = 1 (log scale)");
			h_stack_zoom->GetXaxis()->SetTitle(kinetic_variable);
			h_stack_zoom->GetYaxis()->SetTitle("# events"); 
			gPad->BuildLegend(0.10,0.76,0.4,0.90,"");
			gPad->SetLogy();

			zoom->Modified();
			zoom->Update();	

		}
		else if (k == 2)
		{
			zoom->Close();
		}
	}
	//To save the plots
	string name_png = string(kinetic_variable) + ".png";
	string name_logy_png = string(kinetic_variable) + "_log.png";
	cnv->Print(name_png.c_str(), "png");
	cnv_logy->Print(name_logy_png.c_str(), "png");
	if (kinetic_variable == kinetic_variables[1]) zoom->Print("zoom_mjj.png","png");

	myapp->Run();

	return 0;

}

