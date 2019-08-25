/*
Comparison of SM vs SM + BSM and SM vs SM + inf for cW = 0.05, 0.1, 0.3(HS), 0.4, 1 (with data stored in ntuples)

c++ -o IntVsBsm1 IntVsBsm1.cpp `root-config --glibs --cflags`
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
	int want_log = 0;
	string input_log = "log";
	if (argc == 1) 
	{
		kinetic_variable = "met";
	}
	else 
	{
		kinetic_variable = argv[1];	
	}
	if (argc > 2 && argv[2] == input_log)
	{
		want_log = 1;
	}
		
		
	TApplication* myapp = new TApplication ("myapp", NULL, NULL);
	TCanvas* cnv = new TCanvas("cnv","cnv",0,0,1200,800);
	TCanvas* cnv2 = new TCanvas("cnv2","cnv2",0,0,1000,800);

	cnv->Divide(3,2);
	cnv2->Divide(2,2);

	string cW[] = {"0p05","0p1","0p3","0p4","1"};
	string titles[] = {"cW = 0.05", "cW = 0.1", "cW = 0.3 (high statistics)","cW = 0.4", "cW = 1"};
	//string words[] = {"ntuple_RcW_",".root","SSeu_RcW_bsm_","SSeu_RcW_int_","_nums","_HS.root"};
	string name_histograms[] = {"SM", "BSM", "interference"};
	string kinetic_variables[] = {"met","mjj","mll","ptl1","ptl2"};
	
	vector<TH1F*> histos;
	
	float max_tot[5];
	float maxima[][5] = {{450, 9000, 900, 600, 250},{600, 9000, 900, 600, 250},
		{600, 9000, 1500, 800, 300},{700, 9000, 1300, 1100, 500},
		{1000, 9000, 2200, 1900, 750}}; //maxima in the histograms (raws--->cWs)
	float RMS_array[5] = {87.1868, 952.47, 115.707, 70.3347, 31.2509}; //for every kinetic_variable, same for different cWs
	int Nbins_array[5];

	for (int i = 0; i < 5; i++) 	
	{
		if (kinetic_variable == kinetic_variables[i])
		{
			for (int j = 0; j < 5; j++) 
			{
				max_tot[j] = maxima[j][i];	
				Nbins_array[j] = floor(max_tot[j]/((1./3.)*RMS_array[i]));
			}
			break;
		}
	}	
	
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
		THStack* h_stack_bsm = new THStack("hs_bsm","");
		THStack* h_stack_int = new THStack("hs_int","");
		
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
			
			TH1F* histo = new TH1F ("histo", name_histograms[j].c_str(), Nbins_array[k], 0., max_tot[k]);
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

			values[j].clear();
			weights[j].clear();	
		}
		
		int color_SM  = kGreen - 8 ;
		int color_INT_BSM = kOrange - 4 ;

		histos[0]->SetFillColor(color_SM);
		histos[1]->SetFillColor(color_INT_BSM);
		histos[2]->SetFillColor(color_INT_BSM);

		histos[1]->SetTitle("SM + BSM");
		histos[2]->SetTitle("SM + interference");
		
		h_stack_bsm->Add(histos[0]);
		h_stack_bsm->Add(histos[1]);

		h_stack_int->Add(histos[0]);
		h_stack_int->Add(histos[2]);
	
		if (k < 3) cnv->cd(k+1);
		else cnv2->cd(k-2);
		
		h_stack_bsm->Draw("HIST");
		TText* T_bsm = new TText(); 
		T_bsm->SetTextFont(42); 
		T_bsm->SetTextAlign(21);
		string xlabel = string(kinetic_variable) + string(" (GeV)");
		h_stack_bsm->GetXaxis()->SetTitle(xlabel.c_str());
		h_stack_bsm->GetYaxis()->SetTitle("# events"); 
		gPad->BuildLegend(0.40,0.70,0.90,0.90,"");
		if (want_log == 1) 
		{
			string title = titles[k] + " (logarithmic scale)";
			T_bsm->DrawTextNDC(.5,.95,title.c_str());
			gPad->SetLogy();
		}
		else T_bsm->DrawTextNDC(.5,.95,titles[k].c_str());

		if (k < 3) cnv->cd(k+4);
		else cnv2->cd(k);

		h_stack_int->Draw("HIST");
		TText* T_int = new TText(); 
		T_int->SetTextFont(42); 
		T_int->SetTextAlign(21);
		h_stack_int->GetXaxis()->SetTitle(xlabel.c_str());
		h_stack_int->GetYaxis()->SetTitle("# events"); 
		gPad->BuildLegend(0.40,0.70,0.90,0.90,"");
		if (want_log == 1) 
		{
			string title = titles[k] + " (logarithmic scale)";
			T_int->DrawTextNDC(.5,.95,title.c_str());
			gPad->SetLogy();
		}
		else T_int->DrawTextNDC(.5,.95,titles[k].c_str());

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
			
		histos.clear();		
	
	}
	//To save the plots
	string name1_png; string name2_png;
	if (want_log != 1)
	{
		name1_png = string(kinetic_variable) + "_IntBsm_1.png";
		name2_png = string(kinetic_variable) + "_IntBsm_2.png";
		//cout << name2_png << endl;
	}
	else 
	{
		name1_png = string(kinetic_variable) + "_IntBsm_1_log.png";
		name2_png = string(kinetic_variable) + "_IntBsm_2_log.png";
	}
	
	cnv->Print(name1_png.c_str(), "png");
	cnv2->Print(name2_png.c_str(), "png");

	myapp->Run();

	return 0;

}

