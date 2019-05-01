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
#include "RooExponential.h"
#include "RooGenericPdf.h"

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

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

#include <string>
#include <vector> 


//a = 427956.953782
//b= 6857.34647382
//c= 1143.687277
//n = 1566
//bg1 = 443
//max = 500000
using namespace std;
using namespace RooFit;
using namespace RooStats;

struct channel
{
  double_t a; 
  double_t b; 
  double_t c; 
  double_t n; 
  double_t bg1; 
  double_t bg2; 
}; 

  channel createChannel(double_t a_1, double_t b_1, double_t c_1, double_t n_1, double_t bg1_1, double_t bg2_1)
  {
    channel newChannel; 
    newChannel.a = a_1; 
    newChannel.b = b_1; 
    newChannel.c = c_1; 
    newChannel.n = n_1; 
    newChannel.bg1 = bg1_1; 
    newChannel.bg2 = bg2_1; 
    return newChannel; 
  }

ModelConfig* makeMyModel (const char* name, RooWorkspace& ws, vector<channel> channels, int numChannel ) {
  // derived from data 
  cout<<"nunChannel: "+to_string(numChannel)<<endl; 
   ws.factory("beta[0,5]");
  for(int i=0; i<numChannel; i++)
    {    
      string num = to_string(i+1);
      string bkg1 = "bg_a"+num+"["+to_string(channels[i].bg1)+"]"; 
      cout<<"background created: "+bkg1<<endl; 
      string bkg2 = "bg_b"+num+"["+to_string(channels[i].bg2)+"]"; 
      cout<<"backgrounf2 created: "+bkg2<<endl; 
      string mu = "EXPR::mu"+num+"('"+to_string(channels[i].a)+"*(beta*beta)+"+to_string(channels[i].b)+"*beta+"+to_string(channels[i].c)+"',beta)";
      cout<<"mu: "+mu<<endl; 
      string obs = "Poisson::pdf_a"+num+"(na"+num+"["+to_string(channels[i].n)+"],sum::mu_a"+num+"(mu"+num+",bg_a"+num+",bg_b"+num+"))"; 
      cout<<"obs: "+obs<<endl; 
      string systematics = "Uniform::l_bg_a"+num+"(1)"; 
      cout<<"systematics: "+systematics<<endl; 
      string model = "PROD::model"+num+"(pdf_a"+num+",l_bg_a"+num+")"; 
      cout<<"model created: "+model<<endl; 
       ws.factory(bkg1.c_str());

       ws.factory(bkg2.c_str()); 

       ws.factory(mu.c_str()); 

       ws.factory(obs.c_str());

       ws.factory(systematics.c_str()); 

       ws.factory(model.c_str()); 
       
       
    }
  cout<<"strings created"<<endl; 
  
  string channelString = "index["; 
  string modelString = "SIMUL:jointModel(index,"; 
  string nuisString = ""; 
  
  for(int i = 0; i<numChannel; i++)
    {
      string num = to_string(i+1);
      if(i != numChannel-1)
	{
	  modelString+="channel"+num+"=model"+num+","; 
	  channelString+="channel"+num+","; 
	  nuisString+="bg_a"+num+",bg_b"+num+","; 
        }
      else
	{
	  modelString+="channel"+num+"=model"+num+")";
	  channelString+="channel"+num+"]"; 
	  nuisString += "bg_a"+num+",bg_b"+num; 
	}
    }
  cout<<"nuisance: "+nuisString<<endl; 
  cout<<"channel: "+channelString<<endl; 
  cout<<"model: "+modelString<<endl; 
  ws.factory(channelString.c_str()); 
  ws.factory(modelString.c_str()); 

  // observables
  ws.defineSet("obs","index");

  // parameters of interest
  ws.defineSet("poi","beta");

  // nuisance parameters
 
  cout << "set defined" << "\n"; 

  ws.defineSet("nuis", nuisString.c_str());

  // ws.factory(Prior.c_str()); 
  ws.factory("Uniform::prior(beta)");

  cout<<"prior set"<<"\n"; 


  // model config
  ModelConfig* modelConfig = new ModelConfig(name);
  modelConfig->SetWorkspace(ws);
  modelConfig->SetPdf("jointModel");
  modelConfig->SetPriorPdf("prior");
  modelConfig->SetParametersOfInterest(*(ws.set("poi")));
  modelConfig->SetNuisanceParameters(*(ws.set("nuis")));
  modelConfig->SetObservables(*(ws.set("obs")));
  ws.import(*modelConfig);
  cout<<"modelConfig created"<<"\n"; 
  return modelConfig;
}

double_t betaToLambda(double_t beta)
{
  return sqrt(1/beta); 
}

void StandardBayesianMCMCDemoBetaCombined2(){ 
  double maxPOI=-999; 
  RooWorkspace* ws = new RooWorkspace("ws");
   channel electronBB = createChannel(427956.953782,6857.34647382,1143.687277,1566,443.93,0); 
  channel muonBE = createChannel(536404.023389,10719.7660236,1495.66405188,2894,1046.22,0);
   channel electronBE  = createChannel(343668.348765,2797.6558426,952.653773138, 1550, 542.75,0);
  channel muonBB = createChannel(351848.272356,6703.69863309,994.323027103 ,1352,403.18,0); 
  vector<channel> multipleChannels = {muonBE, muonBB, electronBB, electronBE}; 
  ModelConfig* mc = makeMyModel ("test", *ws, multipleChannels, multipleChannels.size());
  RooDataSet data ("data","",*(mc->GetObservables()));
  //ws->var("na")->setVal(7); 
  data.add( *(mc->GetObservables()));
  ws->import (data); // not really needed for your macro

  RooAbsReal* nll = (mc->GetPdf())->createNLL(data); 
  // RooMinuit m(*nll);
  // m.setVerbose(kTRUE);

  SequentialProposal sp(0.1); 
  MCMCCalculator mcmc(data, *mc); 
  cout<<"mcmc created"<<"\n"; 
  mcmc.SetConfidenceLevel(0.95);
  // mcmc.SetVerbose(true); 
 
  mcmc.SetLeftSideTailFraction(0.); 
  mcmc.SetProposalFunction(sp);
  mcmc.SetNumIters(1000000);         // Metropolis-Hastings algorithm iterations, 1000000
  mcmc.SetNumBurnInSteps(50);       // first N steps to be ignored as burn-in

    RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
  if (maxPOI != -999) 
     firstPOI->setMax(maxPOI);
  cout << "POI set" << "\n";  
  
  MCMCInterval* interval = mcmc.GetInterval();
    if (!interval) { 
     cout << "Error computing Bayesian interval - exit " << endl;
     return;
  }
  cout << "MCMC interval created!" << "\n"; 

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

  double_t lowerLimitLambda = betaToLambda(interval->LowerLimit(*firstPOI)); 
  double_t upperLimitLambda = betaToLambda(interval->UpperLimit(*firstPOI)); 

  cout << "\n95% interval on lambda  is : ["<<
    lowerLimitLambda  << ", "<<
    upperLimitLambda  <<"] "<<endl;

  /* get Bayesian Limit
  double cl95Bayesian = BayesianUpperLimit (data
  cout << "Bayesian UL: " << cl95Bayesian << endl;
  */

  // clean up
  delete ws;
  delete mc;

}
