/*
Comparison of SM vs SM + BSM and SM vs SM + inf for cW = 0.05, 0.1, 0.3(HS), 0.4, 1 (with data stored in ntuples)

c++ -o IntVsBsm2 IntVsBsm2.cpp `root-config --glibs --cflags`
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
	TCanvas* cnv2 = new TCanvas("cnv2","cnv2",0,0,1000,400);
	TCanvas* cnv3 = new TCanvas("cnv3","cnv3",0,0,700,500);

	cnv->Divide(3,1);
	cnv2->Divide(2,1);

	string cW[] = {"0p05","0p1","0p3","0p4","1"};
	string titles[] = {"cW = 0.05", "cW = 0.1", "cW = 0.3 (high statistics)","cW = 0.4", "cW = 1"};
	string name_histograms[] = {"SM", "BSM", "interference"};
	string kinetic_variables[] = {"met","mjj","mll","ptl1","ptl2"};
	
	vector<TH1F*> histos;
	
	float maxima[5] = {1000, 9000, 2200, 1900, 750}; //same for all cWs because the distributions are superimosed
	float max_tot;

	float RMS_array[5] = {87.1868, 952.47, 115.707, 70.3347, 31.2509}; //for every kinetic_variable, same for different cWs
	int Nbins_array[5];

	for (int i = 0; i < 5; i++) 	
	{
		if (kinetic_variable == kinetic_variables[i])
		{
			max_tot = maxima[i];
			for (int j = 0; j < 5; j++) 
			{	
				Nbins_array[j] = floor(max_tot/((1./3.)*RMS_array[i]));
			}
			break;
		}
	}	
	
	THStack* h_stack = new THStack("hs","");
	int colors[5] = {kGreen-8, kOrange-4, kCyan-3, kBlue-6, kRed-3};
	
	for (int k = 0; k < 5; k++) // cW = 0.05, 0.1, 0.3, 0.4, 1.
	{
		string name_files[3];
		if (k != 2) //for cW = 0.3 high statistics (HS)
		{
			name_files[0] = "ntuple_SMlimit_HS.root";
			name_files[1] = "ntuple_RcW_" + cW[k] + ".root";
			name_files[2] = "ntuple_RcW_" + cW[k] + ".root";
		}
		else 
		{
			name_files[0] = "ntuple_SMlimit_HS.root";
			name_files[1] = "ntuple_RcW_" + cW[k] + "_HS.root";
			name_files[2] = "ntuple_RcW_" + cW[k] + "_HS.root";
		}

		string name_ntuples[] = {"SSeu_SMlimit", "SSeu_RcW_bsm_" + cW[k], "SSeu_RcW_int_" + cW[k]};
		string name_global_numbers[] = {"SSeu_SMlimit_nums", "SSeu_RcW_bsm_"+ cW[k] +"_nums", "SSeu_RcW_int_" + cW[k] + "_nums"};

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
		
			TH1F* histo = new TH1F ("histo", name_histograms[j].c_str(), Nbins_array[k], 0., max_tot);			
				
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
			
			if (j == 2) 
			{
				for (int n = 1; n <= histo->GetNbinsX(); n++)
				{
					if (histo->GetBinContent(n) < 0) histo->SetBinContent(n,-histo->GetBinContent(n));
				}
			}

			histo->Scale(normalization);
			histos.push_back(histo);	

			values[j].clear();
			weights[j].clear();	
		}
		
		int color_SM  = kGreen - 8 ;
		int color_INT_BSM = kOrange - 4 ;

		histos[0]->SetLineColor(colors[k]);
		histos[1]->SetLineColor(colors[k]);
		histos[2]->SetLineColor(colors[k]);		

		TH1F* BSM_minus_int = new TH1F(*histos[1]);
		histos[2]->Scale(-1);
		BSM_minus_int->Add(histos[2]);
		BSM_minus_int->SetLineWidth(2);

		h_stack->Add(BSM_minus_int);

		string title = string("BSM - abs(interference) for ") + titles[k];
		BSM_minus_int->SetTitle(title.c_str());
	
		if (k < 3) cnv->cd(k+1);
		else cnv2->cd(k-2);
		
		BSM_minus_int->Draw("HIST");
		string xlabel = string(kinetic_variable) + string(" (GeV)");
		BSM_minus_int->GetXaxis()->SetTitle(xlabel.c_str());
		BSM_minus_int->GetYaxis()->SetTitle("# events"); 
			
		histos.clear();		
	
	}

	cnv3->cd();

	h_stack->Draw("hist nostack");
	gPad->SetLogy();
	TText* T = new TText(); 
	T->SetTextFont(42); 
	T->SetTextAlign(21);
	string xlabel = string(kinetic_variable) + string(" (GeV)");
	h_stack->GetXaxis()->SetTitle(xlabel.c_str());
	h_stack->GetYaxis()->SetTitle("# events"); 
	gPad->BuildLegend(0.40,0.70,0.90,0.90,"");
	T->DrawTextNDC(.5,.95,"BSM - abs(interference) --- Log scale");

	cnv3->Update();
	cnv3->Modified();
	
	//To save the plots
	string name1_png = string(kinetic_variable) + "_IntBsm2_1.png";
	string name2_png = string(kinetic_variable) + "_IntBsm2_2.png";
	string name3_png = string(kinetic_variable) + "_IntBsm2_tot.png";
	
	cnv->Print(name1_png.c_str(), "png");
	cnv2->Print(name2_png.c_str(), "png");
	cnv3->Print(name3_png.c_str(), "png");

	myapp->Run();

	return 0;

}

