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

#include <iterator>
#include <map>

#include <string>
#include <vector> 


using namespace std;
using namespace RooFit;
using namespace RooStats;

//user-defined type for the three types of systematics that exist
enum typeSystematic {SIGNAL, JET, OTHER};

//lowerScale and upperScale allow systematics to be asymmetric
struct SystematicValues
{
  double upperScale;
  double lowerScale; 
  SystematicValues(double lower, double upper):lowerScale(lower), upperScale(upper) {}
  SystematicValues(double scale):SystematicValues(scale, scale) {}
  SystematicValues():SystematicValues(1) {}
  
  
};

//overloads the << operator for the SystematicValues struct
ostream& operator<<(ostream& os, const SystematicValues& sysVals)
  {
    os << sysVals.lowerScale << "/" << sysVals.upperScale;
    return os;
  }

//creates the Channel object with the parameters of an individual channel
struct Channel
{
  double_t a; 
  double_t b; 
  double_t c; 
  double_t n; 
  double_t bg1; 
  double_t bg2;
  //maps the type of the systematic to the vector of the systematic values; map should contain only three pairs (1 for each type of systematic)
  map<typeSystematic, vector<SystematicValues>> systematics; 
}; 


ModelConfig* makeMyModel (const char* name, RooWorkspace& ws, vector<Channel> channels, int numChannel) {
  //This for loop outputs the values of systematics for each channel; not necessary for functonality of code, only useful for debugging purposes
  for(int j = 0; j<numChannel; j++)
    {
      for(map<typeSystematic,vector<SystematicValues>>::iterator it = (channels[j].systematics).begin(); it != (channels[j].systematics).end(); it++)
	{
	  typeSystematic type = it->first; 
	  switch(type)
	    {
	    case SIGNAL:
	      for(int i = 0; i<(it->second).size(); i++)
		{
		  cout << "Signal: " << (it->second)[i] << endl; 
		}
	      break; 
	    case JET:
	      for(int i = 0; i<(it->second).size(); i++)
		{
		  cout << "Jet: " << (it->second)[i] << endl; 
		}
	      break; 
	    case OTHER:
	      for(int i = 0; i<(it->second).size(); i++)
		{
		  cout << "Other: " << (it->second)[i] << endl; 
		}
	      break; 
	    }
	}
    }
  //creates a vector of strings for the systematics
  vector<string> systematicNames = {"lumi", "res", "massScale", "zPeak", "trig", "jets", "xSecOther", "pdf", "id", "stats", "pu"}; 
  //creates the Gaussian distributions for the systematics
  for(int i = 0; i<11; i++)
    {
      string gaussianString = "Gaussian::gauss"+systematicNames[i]+"(z[-5,5], mean[0], sigma[1])"; 
      ws.factory(gaussianString.c_str()); 
    }
  ws.factory("beta[0,5]");
  //creates backgrounds and model parameter mu for each channel
  for(int i=0; i<numChannel; i++)
    {    
      string num = to_string(i+1);
      string bkg1 = "bg_a"+num+"["+to_string(channels[i].bg1)+"]"; 
      cout<<"background created: "+bkg1<<endl; 
      string bkg2 = "bg_b"+num+"["+to_string(channels[i].bg2)+"]"; 
      cout<<"background2 created: "+bkg2<<endl; 
      string mu = "EXPR::mu"+num+"('"+to_string(channels[i].a)+"*(beta*beta)+"+to_string(channels[i].b)+"*beta+"+to_string(channels[i].c)+"',beta)";
      cout<<"mu: "+mu<<endl; 

       ws.factory(bkg1.c_str());

       ws.factory(bkg2.c_str()); 
       
       ws.factory(mu.c_str()); 
    } 
  //incorporates systematics for each channel (THIS PART IS NOT FULLY FUNCTIONAL)
  for(int j = 0; j<numChannel; j++)
    {
      string num = to_string(j+1); //creates a number label for each variable, so that different channels don't share variables with the same name
      //loops over the map of systematic values
      for(map<typeSystematic,vector<SystematicValues>>::iterator it = (channels[j].systematics).begin(); it != (channels[j].systematics).end(); it++)
	{
	  typeSystematic type = it->first;
	  //filters systematics by type
	  switch(type)
	    {
	    case SIGNAL:
	      { 
		//multiply the Gaussian and add 1 for each systematic
		for(int i = 0; i<(it->second).size(); i++)
		  {
		    //multiply the Gaussian distribution of the chosen systematic by the actual value of the systematic stored in the map
		    string systematicGaussian = "PROD::syst"+systematicNames[i]+"(gauss"+systematicNames[i]+","+systematicNames[i]+"["+to_string((it->second)[i].upperScale)+"])";
		    ws.factory(systematicGaussian.c_str());  
		    //add 1 to the Gaussian distribution of the systematic
		    string addSystematics = "EXPR::newSyst"+systematicNames[i]+"('syst"+systematicNames[i]+"+1',syst"+systematicNames[i]+")";
		    ws.factory(addSystematics.c_str()); 
		  }
		//multiply the model parameter mu by every systematic (represented by the Gaussian distribution)
		string multiplySystematics = "PROD::newMu"+num+"(mu"+num+","; 
		for(int k = 0; k<(it->second).size()-1; k++)
		  {
		    multiplySystematics+="newSyst"+systematicNames[k]+","; 
		  }
		multiplySystematics+="newSyst"+systematicNames[systematicNames.size()-1]+")";
		ws.factory(multiplySystematics.c_str()); 
	      }
	      break; 
	      //JET case does not work yet
	    case JET:
	      for(int i = 0; i<(it->second).size(); i++)
		{

		}
	      break; 
	      //OTHER case does not work yet
	    case OTHER:
	      for(int i = 0; i<(it->second).size(); i++)
		{

		}
	      break; 
	    }
	}
    }
  //for each channel, create a Poisson distribution with mu and n
  for(int i=0; i<numChannel; i++)
    {
      string num = to_string(i+1); 
      string obs = "Poisson::pdf_a"+num+"(na"+num+"["+to_string(channels[i].n)+"],sum::mu_a"+num+"(newMu"+num+",bg_a"+num+",bg_b"+num+"))"; 
      cout<<"obs: "+obs<<endl;  
      string model = "PROD::model"+num+"(pdf_a"+num+",Uniform::(1))"; 
      cout<<"model created: "+model<<endl; 

       ws.factory(obs.c_str());

       ws.factory(model.c_str()); 
       
       
    }
  //combine all of the different channels with the SIMUL function 
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


  // create the model config
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

//converts beta to lambda since beta = 1/(lambda)^2
double_t betaToLambda(double_t beta)
{
  return sqrt(1/beta); 
}

void StandardBayesianMCMCDemoBetaCombinedSystematics(){ 
  double maxPOI=-999;
  
  //create a RooStats Workspace
  RooWorkspace* ws = new RooWorkspace("ws");

  //create three vectors of systematics, one for each type of systematic
  vector<SystematicValues> signalSystematics = {{1.025}, {.960}, {.988}, {}, {1.003}, {}, {}, {1.075}, {1.001}, {1.120}, {1.011, 0.982}}; 
  vector<SystematicValues> otherSystematics = {{}, {0.990}, {1.006}, {1.050}, {1.003}, {1.500}, {1.070}, {1.171}, {1.000}, {1.017}, {.998,1.001}}; 
  vector<SystematicValues> jetSystematics = {{}, {}, {}, {}, {}, {1.500}, {}, {}, {}, {}, {}};

  //create a map of the systematics for channel 1 (with each additional channel, you must create another map)
  map<typeSystematic, vector<SystematicValues>> channel1Systematics; 

  //add the vectors of systematics of each type to the map of systematics
  channel1Systematics.insert(pair<typeSystematic,vector<SystematicValues>>(SIGNAL, signalSystematics));  
  channel1Systematics.insert(pair<typeSystematic,vector<SystematicValues>>(OTHER, otherSystematics)); 
  channel1Systematics.insert(pair<typeSystematic,vector<SystematicValues>>(JET, jetSystematics));

  //creates the first channel object
  Channel channelOne = {427956.953782,6857.34647382,1143.687277,1566,443.93,0,channel1Systematics}; 

  //store all channels in a vector
  vector<Channel> multipleChannels= {channelOne};  

  //call the "makeMyModel" function to create the model config
  ModelConfig* mc = makeMyModel ("test", *ws, multipleChannels, multipleChannels.size());

  //the rest of the code has remained mostly the same
  RooDataSet data ("data","",*(mc->GetObservables()));
  data.add( *(mc->GetObservables()));
  ws->import (data);
  RooAbsReal* nll = (mc->GetPdf())->createNLL(data); 
  SequentialProposal sp(0.1); 
  MCMCCalculator mcmc(data, *mc); 
  cout<<"mcmc created"<<"\n"; 
  mcmc.SetConfidenceLevel(0.95);
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
  int iPad=1;
  while( (nuis = (RooRealVar*) it->Next() )){
    c2->cd(iPad++);
    plot.DrawChainScatter(*firstPOI,*nuis);
  }

  // print out the iterval on the first Parameter of Interest
  cout << "\n95% interval on " <<firstPOI->GetName()<<" is : ["<<
    interval->LowerLimit(*firstPOI) << ", "<<
    interval->UpperLimit(*firstPOI) <<"] "<<endl;

  //calculates lambda after beta is calculated
  double_t lowerLimitLambda = betaToLambda(interval->LowerLimit(*firstPOI)); 
  double_t upperLimitLambda = betaToLambda(interval->UpperLimit(*firstPOI)); 

  cout << "\n95% interval on lambda  is : ["<<
    lowerLimitLambda  << ", "<<
    upperLimitLambda  <<"] "<<endl;

  // clean up
  delete ws;
  delete mc;

}
