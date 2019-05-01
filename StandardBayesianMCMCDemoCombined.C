#include "TStopwatch.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooAbsReal.h"

#include "RooStats/RooStatsUtils.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SequentialProposal.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

ModelConfig* makeMyModel (const char* name, RooWorkspace& ws) {
  //create beta
  ws.factory("beta[1,0,20]"); 
  
  // derived from data
  ws.factory("sig1[1]"); // POI

  //multiply beta and signal
  ws.factory("PROD::sigbeta1(sig1,beta)"); 

  // predefined nuisances
  ws.factory("bg_a1[6.5]");
  //ws.factory("bg_a1[0]");

  ws.factory("Poisson::pdf_a1(na1[7],sum::mu_a1(sigbeta1,bg_a1))");

  // nuisance PDFs (systematics)
  // ws.factory("Lognormal::l_bg_a1(bg_a1,nom_bg_a1[5,0,20],sum::kappa_bg_a(1,d_bg_a[0.1]))");
  ws.factory("Uniform::l_bg_a1(1)"); 

  // model
  ws.factory("PROD::model1(pdf_a1,l_bg_a1)");

  ws.factory("sig2[1]"); // POI

  //multiply beta and signal
   ws.factory("PROD::sigbeta2(sig2,beta)"); 

  // predefined nuisances
  ws.factory("bg_a2[1.8]");
  //ws.factory("bg_a2[0]");

  ws.factory("Poisson::pdf_a2(na2[2],sum::mu_a2(sigbeta2,bg_a2))");

  // nuisance PDFs (systematics)
  //ws.factory("Lognormal::l_bg_a2(bg_a2,nom_bg_a2[5,0,20],kappa_bg_a)");
  ws.factory("Uniform::l_bg_a2(1)"); 

  // model
  ws.factory("PROD:model2(pdf_a2,l_bg_a2)");
  cout << "second model created!" << "\n"; 

  //Create discrete observable to label channels
   ws.factory("index[channel1,channel2]"); 

   // ws.factory("PROD::jointModel(model1,model2)"); 
  
  //combine models
   ws.factory("SIMUL:jointModel(index,channel1=model1,channel2=model2)"); 
  cout << "jointModel created!" << "\n"; 

  // observables
  ws.defineSet("obs","index");

  // parameters of interest
  ws.defineSet("poi","beta");

  // nuisance parameters
  ws.defineSet("nuis","bg_a1,bg_a2");
  cout << "nuisance parameters set!" << "\n"; 

  // prior (for Bayesian calculation)
  ws.factory("Uniform::prior(beta)");
  /*
  RooPlot* frame1 = ws.var("beta")->frame(); 
  ws.pdf("model1")->plotOn(frame1); 
  frame1->Draw(); 
  */ 
   // model config
  ModelConfig* modelConfig = new ModelConfig(name);
  modelConfig->SetWorkspace(ws);
  modelConfig->SetPdf("jointModel");
  modelConfig->SetPriorPdf("prior");
  modelConfig->SetParametersOfInterest(*(ws.set("poi")));
  modelConfig->SetNuisanceParameters(*(ws.set("nuis")));
  modelConfig->SetObservables(*(ws.set("obs")));
  ws.import(*modelConfig);
  return modelConfig;
}

void StandardBayesianMCMCDemoCombined(){
  double maxPOI=-999; 
  RooWorkspace* ws = new RooWorkspace("ws");
  ModelConfig* mc = makeMyModel ("test", *ws);
  cout << "Model created!" << "\n"; 
  RooDataSet data ("data","",*(mc->GetObservables()));
  data.add( *(mc->GetObservables()));
  ws->import (data); // not really needed for your macro

  RooPlot* frame = ws->var("beta")->frame(); 
  ws->pdf("jointModel")->plotOn(frame, ProjWData(data)); 
  frame->Draw();

  SequentialProposal sp(0.1); 

  MCMCCalculator mcmc(data, *mc);
  cout << "MCMC Calculator created!" << "\n"; 
  mcmc.SetConfidenceLevel(0.95);
 
  mcmc.SetLeftSideTailFraction(0.); 
  mcmc.SetProposalFunction(sp);
  mcmc.SetNumIters(1000000);         // Metropolis-Hastings algorithm iterations
  mcmc.SetNumBurnInSteps(50);       // first N steps to be ignored as burn-in
  cout << "MCMC Things set!" << "\n"; 
    RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
  if (maxPOI != -999) 
     firstPOI->setMax(maxPOI);
  cout << "POI set!" << "\n"; 

  MCMCInterval* interval = mcmc.GetInterval();
  cout << "MCMC interval created!" << "\n"; 

  // make a plot
  //TCanvas* c1 =
  new TCanvas("IntervalPlot");
  MCMCIntervalPlot plot(*interval);
  plot.Draw();
  cout << "Plot drawn!" << "\n"; 

  TCanvas* c2 = new TCanvas("extraPlots");
  const RooArgSet* list = mc->GetNuisanceParameters();
  if(list->getSize()>1){
    double n = list->getSize();
    int ny = TMath::CeilNint( sqrt(n) );
    int nx = TMath::CeilNint(double(n)/ny);
    c2->Divide( nx,ny);
  }

  // draw a scatter plot of chain results for poi vs each nuisance parameters
  TIterator* it = mc->GetNuisanceParameters()->createIterator();
  RooRealVar* nuis = NULL;
  int iPad=1; // iPad, that's funny
  while( (nuis = (RooRealVar*) it->Next() )){
    c2->cd(iPad++);
    plot.DrawChainScatter(*firstPOI,*nuis);
  }


  // print out the iterval on the first Parameter of Interest
  cout << "\n95% interval on " <<firstPOI->GetName()<<" is : ["<<
    interval->LowerLimit(*firstPOI) << ", "<<
    interval->UpperLimit(*firstPOI) <<"] "<<endl;


  // clean up
  delete ws;
  delete mc;

}

 

 
 
