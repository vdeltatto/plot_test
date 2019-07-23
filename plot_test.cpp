/* 
   root[] .L plot_test.cpp
   root[] plot_test("kinematic_variable")
*/

#include <TFile.h>
#include <TNtuple.h>
#include <TTreeReader.h>
#include <TCanvas.h>
#include <TGedEditor.h>

using namespace std;


void plot_test (const char* kinematic_variable = "met") //other variables: mjj, mll, ptl1, ptl2, w
{
	
	TFile* myfile = new TFile ("ntuple_SMlimit.root") ;
	TNtuple* nt = (TNtuple*) myfile->Get("SSeu_SMlimit") ;

	TCanvas* cnv = new TCanvas(kinematic_variable,kinematic_variable,0,0,1200,400) ;
	 
	cnv->Divide (3,1);
	cnv->cd (1) ;  
	nt->Draw(kinematic_variable) ;
	nt->GetHistogram()->SetTitle("SM simulation");
	nt->GetHistogram()->GetXaxis()->SetTitle(kinematic_variable);
	//nt->GetHistogram()->GetYaxis()->SetTitle("# events");
	nt->GetHistogram()->SetFillColor(kGreen);
	cnv->Modified () ;
	cnv->Update () ; 
	

	myfile->Close() ;
	myfile = new TFile ("ntuple_RcW_0p3.root") ;
	nt = (TNtuple*) myfile->Get("SSeu_RcW_bsm_0p3") ;

	cnv->cd (2) ;  
	nt->Draw(kinematic_variable) ;
	nt->GetHistogram()->SetTitle("BSM (quadratic term)");
	nt->GetHistogram()->GetXaxis()->SetTitle(kinematic_variable);
	//nt->GetHistogram()->GetYaxis()->SetTitle("# events");
	nt->GetHistogram()->SetFillColor(kBlue);
	cnv->Modified () ;
	cnv->Update () ; 
	

	nt = (TNtuple*) myfile->Get("SSeu_RcW_int_0p3") ;

	cnv->cd (3) ;  
	nt->Draw(kinematic_variable) ;
	nt->GetHistogram()->SetTitle("BSM (interference term)");
	nt->GetHistogram()->GetXaxis()->SetTitle(kinematic_variable);
	//nt->GetHistogram()->GetYaxis()->SetTitle("# events");
	nt->GetHistogram()->SetFillColor(kRed);
	cnv->Modified () ;
	cnv->Update () ; 
	

}
