/*
c++ -o createFiles createFiles.cpp `root-config --glibs --cflags` 
*/

#include <iostream>
#include <string>
#include <iomanip>

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
	
	TH1::SetDefaultSumw2();
	string wilson_coefficients[] = {"0p05","0p1","0p3","0p4","1"};

	const char* cW;

	if (argc == 1) 
	{
		cW = "0p3";
	}
	else 
	{
		cW = argv[1];
	}

	string name_f_met = "histo_" + string(cW) + "_met.root";
	string name_f_mjj = "histo_" + string(cW) + "_mjj.root";
	string name_f_mll = "histo_" + string(cW) + "_mll.root";
	string name_f_ptl1 = "histo_" + string(cW) + "_ptl1.root";
	string name_f_ptl2 = "histo_" + string(cW) + "_ptl2.root";

	vector<TFile*> files;
	TFile* f_met = new TFile (name_f_met.c_str(),"recreate");
	TFile* f_mjj = new TFile (name_f_mjj.c_str(),"recreate");
	TFile* f_mll = new TFile (name_f_mll.c_str(),"recreate");
	TFile* f_ptl1 = new TFile (name_f_ptl1.c_str(),"recreate");
	TFile* f_ptl2 = new TFile (name_f_ptl2.c_str(),"recreate");

	files.push_back(f_met);
	files.push_back(f_mjj);
	files.push_back(f_mll);
	files.push_back(f_ptl1);
	files.push_back(f_ptl2);
	
	const char* kinetic_variable; //possible variables: met, mjj, mll, ptl1, ptl2

	
	string name_histograms[] = {"histo_sm", "histo_linear", "histo_quadratic"};
	string kinetic_variables[] = {"met","mjj","mll","ptl1","ptl2"};
	
	string name_files[3];
	if (cW != wilson_coefficients[2]) //for cW = 0.3 high statistics (HS)
		{
			name_files[0] = "ntuple_SMlimit_HS.root";
			name_files[1] = "ntuple_RcW_" + string(cW) + ".root";
			name_files[2] = "ntuple_RcW_" + string(cW) + ".root";
		}
		else 
		{
			name_files[0] = "ntuple_SMlimit_HS.root";
			name_files[1] = "ntuple_RcW_" + string(cW) + "_HS.root";
			name_files[2] = "ntuple_RcW_" + string(cW) + "_HS.root";
		}
	string name_ntuples[] = {"SSeu_SMlimit", "SSeu_RcW_int_" + string(cW), "SSeu_RcW_bsm_" + string(cW)};
	string name_global_numbers[] = {"SSeu_SMlimit_nums", "SSeu_RcW_int_"+string(cW)+"_nums", "SSeu_RcW_bsm_" + string(cW) + "_nums"};

	float max[5];
	float maxima[][5] = {{1600, 7500, 1900, 1200, 600},{1600, 7500, 1900, 1200, 600}, 
		{1600, 7500, 1900, 1200, 600},{1600, 7500, 1900, 1200, 600},
		{1600, 7500, 1900, 1200, 600}};

	float RMS_array[5] = {87.1868, 952.47, 115.707, 70.3347, 31.2509}; //for every kinetic variable (same for different cWs)

	//every raw for every cW
	//int Nbins_array[][5] = {{25, 26, 34, 36, 35},{31, 26, 34, 36, 35},{28, 26, 37, 35, 40},{33, 26, 50, 52, 54},{38, 26, 54, 48, 58}};
	/*int Nbins_array[5][5];
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			Nbins_array[i][j] = floor(maxima[i][j]/(0.25*RMS_array[j]));	
			//cout << Nbins_array[i][j] << endl;
		}
	}*/
	int Nbins[5];

	for (int i = 0; i < 5; i++) 	
	{
		if (cW == wilson_coefficients[i])
		{
			for (int j = 0; j < 5; j++)
			{
				max[j] = maxima[i][j];
				Nbins[j] = floor(max[j]/((1./3.)*RMS_array[j]));	
			}
			break;
		}
	}	

	vector<float> values; 
	vector<float> weights;
	vector<TH1F*> histos;
	cout << setprecision(7) << fixed;
		
	for (int k = 0; k < 5; k++) 
	{
		cout << kinetic_variables[k].c_str() << "---------------------------------" << endl;
		kinetic_variable = kinetic_variables[k].c_str();
		for (int j = 0; j < 3; j++) // j = 0,1,2: SM simulation, BSM (interference term), BSM (quadratic term)
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
			
			/*min = *min_element(values.begin(), values.end());
			max = *max_element(values.begin(), values.end());*/
	
			TH1F* histo = new TH1F (name_histograms[j].c_str(), kinetic_variable, Nbins[k], 0., max[k]);		
				
			TH1F* global_numbers = (TH1F*) myfile->Get(name_global_numbers[j].c_str()) ;
			float cross_section = global_numbers->GetBinContent(1);
			float sum_weights_total = global_numbers->GetBinContent(2);
			float sum_weights_selected = global_numbers->GetBinContent(3);
			float luminosity = 100;
			float normalization = cross_section*luminosity/sum_weights_total;

			files[k]->cd();

			for (int i = 0; i < values.size(); i++) 
			{
				histo->Fill(values[i],weights[i]);
			}
	
			histo->Scale(normalization);
			histo->SetBinContent(histo->GetNbinsX(), histo->GetBinContent(histo->GetNbinsX()) + histo->GetBinContent(histo->GetNbinsX() + 1));
			histo->SetBinContent(histo->GetNbinsX() + 1, 0.);

			/*for (int i = 0; i < Nbins[k] + 2 ; i++) {
				if (histo->GetBinContent(i) != 0)
					cout << histo->GetBinContent(i) << endl;
				else
					cout << "XXXXXXX" << endl;
			}
			cout << endl;*/

			histo->Write();
			cout << histo->Integral() << "\t";
			files[k]->Write();

			values.clear();
			weights.clear();	
		}
		files[k]->Close();
		cout << endl << endl;
	}

	return 0 ;
}

