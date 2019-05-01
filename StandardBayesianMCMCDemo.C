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
  // derived from data
  ws.factory("sig[2,0,20]"); // POI

  // predefined nuisances
  ws.factory("bg_a[2]");

  ws.factory("Poisson::pdf_a(na[2],sum::mu_a(sig,bg_a))");

  // nuisance PDFs (systematics)
  // ws.factory("Gaussian:l_bg_a(bg_a,bg_mean[0.5],sigma[0.25])"); 
  ws.factory("Uniform::l_bg_a(1)"); 

  // model 
  ws.factory("PROD::model(pdf_a,l_bg_a)");

  // observables
  ws.defineSet("obs","na");

  // parameters of interest
  ws.defineSet("poi","sig");

  // nuisance parameters
  ws.defineSet("nuis","bg_a");

  // prior (for Bayesian calculation)
  ws.factory("Uniform::prior(sig)");


  // model config
  ModelConfig* modelConfig = new ModelConfig(name);
  modelConfig->SetWorkspace(ws);
  modelConfig->SetPdf("model");
  modelConfig->SetPriorPdf("prior");
  modelConfig->SetParametersOfInterest(*(ws.set("poi")));
  modelConfig->SetNuisanceParameters(*(ws.set("nuis")));
  modelConfig->SetObservables(*(ws.set("obs")));
  ws.import(*modelConfig);
  return modelConfig;
}

double BayesianUpperLimit (RooAbsData& data, ModelConfig& modelConfig, double CL = 0.95) {
  BayesianCalculator calculator (data, modelConfig);
  calculator.SetConfidenceLevel(CL);
  calculator.SetLeftSideTailFraction (0.); // UL
  SimpleInterval* interval = calculator.GetInterval();
  calculator.SetScanOfPosterior(40);
  RooPlot * plot = calculator.GetPosteriorPlot();
  plot->Draw(); 
  return interval->UpperLimit();
}

void StandardBayesianMCMCDemo(){
  double maxPOI=-999; 
  RooWorkspace* ws = new RooWorkspace("ws");
  ModelConfig* mc = makeMyModel ("test", *ws);
  RooDataSet data ("data","",*(mc->GetObservables()));
  //ws->var("na")->setVal(7); 
  data.add( *(mc->GetObservables()));
  ws->import (data); // not really needed for your macro

  SequentialProposal sp(0.1); 

  MCMCCalculator mcmc(data, *mc); 
  mcmc.SetConfidenceLevel(0.95);
 
  mcmc.SetLeftSideTailFraction(0.); 
  mcmc.SetProposalFunction(sp);
  mcmc.SetNumIters(1000000);         // Metropolis-Hastings algorithm iterations
  mcmc.SetNumBurnInSteps(50);       // first N steps to be ignored as burn-in

    RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
  if (maxPOI != -999) 
     firstPOI->setMax(maxPOI);

  MCMCInterval* interval = mcmc.GetInterval();

  // make a plot
  //TCanvas* c1 =
  new TCanvas("IntervalPlot");
  MCMCIntervalPlot plot(*interval);
  plot.Draw();

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

  /* get Bayesian Limit
  double cl95Bayesian = BayesianUpperLimit (data
  cout << "Bayesian UL: " << cl95Bayesian << endl;
  */

  // clean up
  delete ws;
  delete mc;

}
