import FWCore.ParameterSet.Config as cms

TTHMultiLepFilter = cms.EDFilter("TTHMultiLepFilter",
                                 genParticles = cms.InputTag("genParticles"),
                                 genJets = cms.InputTag("ak4GenJetsCustom"),
                                 )
