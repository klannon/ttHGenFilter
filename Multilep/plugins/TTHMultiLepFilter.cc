// -*- C++ -*-
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"


//
// class declaration
//

class TTHMultiLepFilter : public edm::stream::EDFilter<> {
public:
  explicit TTHMultiLepFilter(const edm::ParameterSet&);
  ~TTHMultiLepFilter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginStream(edm::StreamID) override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;
  
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  const edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
  
  
  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TTHMultiLepFilter::TTHMultiLepFilter(const edm::ParameterSet& iConfig):
genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))) {

}


TTHMultiLepFilter::~TTHMultiLepFilter() {

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool TTHMultiLepFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   //get GenParticleCollection
   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genParticlesToken_, genParticles);

   //get GenJetCollection
   edm::Handle<reco::GenJetCollection> genJets;
   iEvent.getByToken(genJetsToken_, genJets);

   //First we need to select the gen leptons
   reco::GenParticleCollection genLeptons;
   for (auto genp : *genParticles) {
     if (genp.status() != 1) continue;
     if (abs(genp.pdgId()) != 11 && abs(genp.pdgId()) != 13) continue;
     if (genp.pt() < 9 || fabs(genp.eta()) > 2.8) continue;

     //If we get here, we have a reasonable lepton.  Now we need to check on its isolation properties

     //Let's start with the jet pt ratio.  Requires looping over all the jets to find a match
     double minDR = 9e20;
     double jetRatio = 1.5;
     for (auto jet: *genJets) {
       if (jet.pt() < 5) continue;
       if (fabs(jet.eta()) > 3) continue;

       double dR = reco::deltaR(genp,jet);
       if (dR < minDR) {
         minDR = dR;
         if (dR < 0.5) jetRatio = std::min(genp.pt()/jet.pt(),1.5f);
       }
     }

     if (jetRatio < 0.6) continue;

     //Final check, miniIso
     double miniIsoR = 10.0/std::min(std::max(float(genp.pt()), float(50.)),float(200.));
     double miniIsoE = 0.;
     for (auto g: *genParticles) {
       //Only Final State Particles
       if (g.status() != 1) continue;

       //Skip neutrinos!
       if (abs(g.pdgId()) == 12 ||
           abs(g.pdgId()) == 14 ||
           abs(g.pdgId()) == 16) continue;
              
       double dR = reco::deltaR(genp,g);
       if (dR < miniIsoR && dR > 0.01 && g.pt() > 0.5) {
         miniIsoE += g.pt();
       }
     }
     
     double miniIso = miniIsoE/genp.pt();

     if (miniIso > 0.3) continue;

     //OK, if we make it to HERE then we've got a lepton that's worth counting
     genLeptons.push_back(genp);

   }


   //Sort by pt
   struct {
     bool operator()(const reco::GenParticle &a, const reco::GenParticle &b) {
       return a.pt() > b.pt();
     }
   } ByPt;
   std::sort(genLeptons.begin(),genLeptons.end(),ByPt);

   //Next step: Check whether our leptons look good
   if (genLeptons.size() < 2) {
     return false;
   } else if (genLeptons.size() == 2) {

     double pt1Cut = (abs(genLeptons[0].pdgId()) == 11) ? 24. : 19.;
     double pt2Cut = (abs(genLeptons[1].pdgId()) == 11) ? 14. :  9.;
     
     if (genLeptons[0].pt() < pt1Cut || genLeptons[1].pt() < pt2Cut) return false;
     //This implements the same-sign requirement
     if (genLeptons[0].pdgId()*genLeptons[1].pdgId() < 0) return false;

   } else { //Now we have 3 or more leptons

     if (genLeptons[1].pt() < 19.) return false;

   }

   // Finally, if we're here, we just need to check on the jets
   int nCleanedJets = 0;
   for (auto jet: *genJets) {
     if (jet.pt() < 5) continue;
     if (fabs(jet.eta()) > 3) continue;

     double minDR = 9e20;
     for (auto lep: genLeptons) {
       double dR = reco::deltaR(jet,lep);
       if (dR < minDR) {
         minDR = dR;
       }
     }

     if (minDR > 0.5 && jet.pt() > 20. && fabs(jet.eta()) < 3.)
       ++nCleanedJets;
   }

   //To get to here, all events are either good 2lss or 3l+ events.
   if (genLeptons.size() == 2) {
     return (nCleanedJets >= 4);
   } else {
     return (nCleanedJets >=2);
   }

}



// ------------ method called once each stream before processing any runs, lumis or events  ------------
void TTHMultiLepFilter::beginStream(edm::StreamID) {

}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void TTHMultiLepFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void TTHMultiLepFilter::beginRun(edm::Run const&, edm::EventSetup const&) {

}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void TTHMultiLepFilter::endRun(edm::Run const&, edm::EventSetup const&) {

}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void TTHMultiLepFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {

}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void TTHMultiLepFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {

}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTHMultiLepFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TTHMultiLepFilter);
