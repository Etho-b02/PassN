//
// analyzer module producing DQM histograms for Stm subsystem
//
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"

#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMFragmentSummary.hh"
#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

#include "TH1F.h"

#include <string>
#include <vector>

namespace mu2e
{

class DqmStm : public art::EDAnalyzer
{
 public:
  struct Config
  {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    //For now using explciit tags
    fhicl::Atom<art::InputTag> rawHPGeTag{Name("rawHPGeTag"), Comment("StmWaveform rawHPGe Collection"), art::InputTag()};
    fhicl::Atom<art::InputTag> rawLaBrTag{Name("rawLaBrTag"), Comment("StmWaveform rawLaBr Collection"), art::InputTag()};
    fhicl::Atom<art::InputTag> zsHPGeTag{Name("zsHPGeTag"), Comment("StmWaveform zsHPGe Collection"), art::InputTag()};
    fhicl::Atom<art::InputTag> zsLaBrTag{Name("zsLaBrTag"), Comment("StmWaveform zsLaBr Collection"), art::InputTag()};
    fhicl::Atom<art::InputTag> phHPGeTag{Name("phHPGeTag"), Comment("StmDigi phHPge Collection"), art::InputTag()};
    fhicl::Atom<art::InputTag> phLaBrTag{Name("phLaBrTag"), Comment("StmDigi phLaBr Collection"), art::InputTag()};
    fhicl::Atom<art::InputTag> HPGefragSummaryTag{Name("HPGefragSummaryTag"), Comment("STM HPGe fragment summary collection"), art::InputTag()};
    fhicl::Atom<art::InputTag> LaBrfragSummaryTag{Name("LaBrfragSummaryTag"), Comment("STM LaBr fragment summary collection"), art::InputTag()};
  };
  typedef art::EDAnalyzer::Table<Config> Parameters;

  explicit DqmStm(const Parameters& conf);
  virtual ~DqmStm() {}

  virtual void beginJob();
  virtual void endJob(){};
  virtual void analyze(const art::Event& e);

private:
  Config _conf;

  //Version number tracked
  TH1D* _hVer;
  // Container and Inner Frag Counts
  TH1D* _hNContainerFragsHPGe;
  TH1D* _hNTotalInnerFragsHPGe;
  TH1D* _hNContainerFragsLaBr;
  TH1D* _hNTotalInnerFragsLaBr;

  // Goal is to show how detector behavior changes over time
  // Our focus is on ph Digis and STMFragSummaries
  //Raw
  TH1D* _hNRAWHPGe; //Raw Waveform size
  TH1D* _hNRAWLaBr;
  TH1D* _hNGoodRawFragsHPGe;
  TH1D* _hNGoodRawFragsLaBr;
  TH1D* _hNZeroRawFragsHPGe;
  TH1D* _hNZeroRawFragsLaBr;
  TH1D* _hNEmptyRawFragsHPGe;
  TH1D* _hNEmptyRawFragsLaBr;
  TH1D* _hRawHPGeADC;
  TH1D* _hRawLaBrADC;
  // We want max adc values per wavelength
  TH1D* _hMaxRawADCHPGe;
  TH1D* _hMaxRawADCLaBr;

  //ZS
  TH1D* _hNZSHPGe; //ZS Waveform size
  TH1D* _hNZSLaBr;
  TH1D* _hNGoodZSFragsHPGe;
  TH1D* _hNGoodZSFragsLaBr;
  TH1D* _hNZeroZSFragsHPGe;
  TH1D* _hNZeroZSFragsLaBr;
  TH1D* _hNEmptyZSFragsHPGe;
  TH1D* _hNEmptyZSFragsLaBr;
  TH1D* _hZSHPGeADC;
  TH1D* _hZSLaBrADC;

  // PH digis
  // TH1D* _hNPHHPGeInnerFrags;
  // TH1D* _hNPHLaBrInnerFrags;
  TH1D* _hNPHHPGe; //Number of pulse heights for HPGe
  TH1D* _hPulseHeightHPGe; //Pulse Height of HPGe
  TH1D* _hNPHLaBr; //Number of pulse heights for LaBr
  TH1D* _hPulseHeightLaBr; //Pulse Height of LaBr

  TH1D* _hNGoodPHFragsHPGe;
  TH1D* _hNGoodPHFragsLaBr;
  TH1D* _hNZeroPHFragsHPGe;
  TH1D* _hNZeroPHFragsLaBr;
  TH1D* _hNEmptyPHFragsHPGe;
  TH1D* _hNEmptyPHFragsLaBr;

};
/*****************************************/
  // Configrues fcl to local variables
DqmStm::DqmStm(const Parameters& conf) : art::EDAnalyzer(conf), _conf(conf()) {
  mayConsume<STMPHDigiCollection>(_conf.phHPGeTag()); //reads PH from the art file given in fcl
  mayConsume<STMPHDigiCollection>(_conf.phLaBrTag());
  mayConsume<STMWaveformDigiCollection>(_conf.rawHPGeTag());
  mayConsume<STMWaveformDigiCollection>(_conf.rawLaBrTag());
  mayConsume<STMWaveformDigiCollection>(_conf.zsHPGeTag());
  mayConsume<STMWaveformDigiCollection>(_conf.zsLaBrTag());
  mayConsume<STMFragmentSummaryCollection>(_conf.HPGefragSummaryTag());
  mayConsume<STMFragmentSummaryCollection>(_conf.LaBrfragSummaryTag());
}
/*****************************************/
  // makes the hist -> working with _hNPH and _hPulseHeights
void DqmStm::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;

  // Histogram update version
  _hVer = tfs->make<TH1D>("Ver",
    "Version Number",
    101, -0.5, 100.00);

  // Need to find a way to run these if we are only saving ph digis and frag summaries
  if (!_conf.rawHPGeTag().empty()){
    _hNRAWHPGe = tfs->make<TH1D>("NRawWaveformDigisHPGe",
      "N RAW HPGe WaveformDigis;Number of Raw WF Digis;Entries",
      250, -0.5, 500.5);

    _hRawHPGeADC = tfs->make<TH1D>("RawADCHPGe",
      "Raw HPGe ADC in Waveform;ADC;Samples",
      500,-500.0,500.0);

    _hMaxRawADCHPGe = tfs->make<TH1D>("MaxADCValuesHPGe",
      "Max ADC per Raw HPGe Waveform; Max ADC; Waveform Count",
      200, 0, 500);
  }

  if (!_conf.rawLaBrTag().empty()){
    _hNRAWLaBr = tfs->make<TH1D>("NRawWaveformDigisLaBr",
      "N RAW LaBr WaveformDigis;Number of Raw WF Digis:Entries",
      250, -0.5, 500.5);

    _hRawLaBrADC = tfs->make<TH1D>("RawADCLaBr",
      "Raw LaBr ADC in Waveform;ADC;Samples",
      500,-500.0,500.0);

    _hMaxRawADCLaBr = tfs->make<TH1D>("MaxADCValuesLaBr",
      "Max ADC per Raw LaBr Waveform; Max ADC; Waveform Count",
      200, 0, 500);
  }

  // Can be potentially used to match up wit ph?
  if (!_conf.zsHPGeTag().empty()){
    _hNZSHPGe = tfs->make<TH1D>("NZSWaveformDigisHPGe",
      "N ZS HPGe WaveformDigis;Number of ZS WF Digis;Entries",
      250, -0.5, 500.5);

    _hZSHPGeADC = tfs->make<TH1D>("ZSADCHPGe",
      "ZS HPGe ADC in Waveform;ADC;Samples",
      500,-500.0,500.0);
  }

  if (!_conf.zsLaBrTag().empty()){
    _hNZSLaBr = tfs->make<TH1D>("NZSWaveformDigisLaBr",
      "N ZS LaBr WaveformDigis;Number of ZS WF Digis;Entries",
      250, -0.5, 500.5);

    _hZSLaBrADC = tfs->make<TH1D>("ZSADCLaBr",
      "ZS LaBr ADC in Waveform;ADC;Samples",
      500,-500.0,500.0);
  }

  // Main DQM focus
  if (!_conf.phHPGeTag().empty()){//Fix this one first, configure the other one later
    _hNPHHPGe = tfs->make<TH1D>("NPHDigisHPGe",
      "N PH HPGe Digis;Number of Pulse Height Digis:Entries",
      20, -0.5, 30.5);

    _hPulseHeightHPGe = tfs->make<TH1D>("PulseHeightHPGe",
      "STM HPGe Pulse Height;Pulse Height;Entries",
      600, -300, 300);
  }

  if (!_conf.phLaBrTag().empty()){
    _hNPHLaBr = tfs->make<TH1D>("NPHDigisLaBr",
      "N PH LaBr Digis;Number of Pulse Height Digis;Entries",
      20, -0.5, 30.5);

      _hPulseHeightLaBr = tfs->make<TH1D>("PulseHeightLaBr",
      "STM LaBr Pulse Height;Pulse Height;Entries",
      600,-300,300);
  }

  if (!_conf.HPGefragSummaryTag().empty()){
    _hNContainerFragsHPGe = tfs->make<TH1D>("NContainerFragsHPGe",
      "HPGe Container Fragments;Number of HPGe Container Frags per art Event;art Events",
      4,-0.5,4);

    _hNTotalInnerFragsHPGe = tfs->make<TH1D>("NTotalInnerFragsHPGe",
      "HPGe Inner Fragments;Number of HPGe Inner Frags per art Event;art Events",
      200,-0.5,2000.5);

    _hNGoodRawFragsHPGe = tfs->make<TH1D>("NGoodRawFragsHPGe",
      "Good Raw Fragments HPGe;Good Raw Fragments per art Event;art Events",
      200,-10.5, 250.5);

    _hNGoodZSFragsHPGe = tfs->make<TH1D>("NGoodZSFragsHPGe",
      "Good ZS Fragments HPGe;Good ZS Fragments per art Event;art Events",
      200,-10.5, 250.5);

    _hNGoodPHFragsHPGe = tfs->make<TH1D>("NGoodPHFragsHPGe",
      "Good PH Fragments HPGe;Good PH Fragments per art Event;art Events",
      200,-10.5,250.5);

    _hNZeroRawFragsHPGe = tfs->make<TH1D>("NZeroRawFragsHPGe",
      "Zero Filled Raw Fragments HPGe;Zero Filled Raw Fragments per art Event;art Events",
      200,-10.5,250.5);

    _hNZeroZSFragsHPGe = tfs->make<TH1D>("NZeroZSFragsHPGe",
      "Zero Filled ZS Fragments HPGe;Zero Filled ZS Fragments per art Event;art Events",
      200,-10.5,250.5);

    _hNZeroPHFragsHPGe = tfs->make<TH1D>("NZeroPHFragsHPGe",
      "Zero Filled PH Fragments HPGe;Zero Filled PH Fragments per art Event;art Events",
      200,-10.5,250.5);

    _hNEmptyRawFragsHPGe = tfs->make<TH1D>("NEmptyRawFragsHPGe",
      "Empty Raw Fragments HPGe;Empty Raw Fragments per art Event;art Events",
      200,-10.5,250.5);

    _hNEmptyZSFragsHPGe = tfs->make<TH1D>("NEmptyZSFragsHPGe",
      "Empty ZS Fragments HPGe;Empty ZS Fragments per art Event;art Events",
      200,-10.5,250.5);

    _hNEmptyPHFragsHPGe = tfs->make<TH1D>("NEmptyPHFragsHPGe",
      "Empty PH Fragments HPGe;Empty PH Fragments per art Event;art Events",
      200,-10.5,250.5);
  }

  if (!_conf.LaBrfragSummaryTag().empty()){
    _hNContainerFragsLaBr = tfs->make<TH1D>("NContainerFragsLaBr",
      "LaBr Container Fragments;Number of LaBr Container Frags per art Event;art Events",
      4,-0.5,4);

    _hNTotalInnerFragsLaBr = tfs->make<TH1D>("NTotalInnerFragsLaBr",
      "LaBr Inner Fragments;Number of LaBr Inner Frgs per art Event ;Events",
      200,-0.5,2000.5);

    _hNGoodRawFragsLaBr = tfs->make<TH1D>("NGoodRawFragsLaBr",
      "Good Raw Fragments LaBr;Good Raw Fragments per art Event;art Events",
      200, -10.5, 250.5);

    _hNGoodZSFragsLaBr = tfs->make<TH1D>("NGoodZSFragsLaBr",
      "Good ZS Fragments LaBr;Good ZS Fragments per art Event;art Events",
      200, -10.5, 250.5);

    _hNGoodPHFragsLaBr = tfs->make<TH1D>("NGoodPHFragsLaBr",
      "Good PH Fragments LaBr;Good PH Fragments per art Event;art Events",
      200, -10.5, 250.5);

    _hNZeroRawFragsLaBr = tfs->make<TH1D>("NZeroRawFragsLaBr",
      "Zero Filled Raw Fragments LaBr;Zero Filled Raw Fragments per art Event;art Events",
      200, -10.5, 250.5);

    _hNZeroZSFragsLaBr = tfs->make<TH1D>("NZeroZSFragsLaBr",
      "Zero Filled ZS Fragments LaBr;Zero Filled ZS Fragments per art Event;art Events",
      200, -10.5, 250.5);

    _hNZeroPHFragsLaBr = tfs->make<TH1D>("NZeroPHFragsLaBr",
      "Zero Filled PH Fragments LaBr;Zero Filled PH Fragments per art Event;art Events",
      200, -10.5, 250.5);

    _hNEmptyRawFragsLaBr = tfs->make<TH1D>("NEmptyRawFragsLaBr",
      "Empty Raw Fragments LaBr;Empty Raw Fragments per art Event;art Events",
      200, -10.5, 250.5);

    _hNEmptyZSFragsLaBr = tfs->make<TH1D>("NEmptyZSFragsLaBr",
      "Empty ZS Fragments LaBr;Empty ZS Fragments per art Event;art Events",
      200, -10.5, 250.5);

    _hNEmptyPHFragsLaBr = tfs->make<TH1D>("NEmptyPHFragsLaBr",
      "Empty PH Fragments LaBr;Empty PH Fragments per art Event;art Events",
      200, -10.5, 250.5);
  }
}

/************************************/
  // Filling the histograms
void DqmStm::analyze(const art::Event& event) {
  _hVer->Fill(0.0);

  if (!_conf.rawHPGeTag().empty()) {
    auto rawDigiHandle = event.getValidHandle<STMWaveformDigiCollection>(_conf.rawHPGeTag());
    const auto& rawdigis = *rawDigiHandle;
    _hNRAWHPGe->Fill(rawdigis.size());

    for (const auto& dg : rawdigis){ //Goes into inner frag
      auto const& adcs = dg.adcs();
      for (const auto& adc: adcs) {
        _hRawHPGeADC->Fill(adc);
      }
      if (adcs.empty()) {continue;}

      auto maxADCH = std::max_element(adcs.begin(),adcs.end());
      _hMaxRawADCHPGe ->Fill(*maxADCH);
  }
}

  if (!_conf.rawLaBrTag().empty()) {
    auto rawDigiHandle = event.getValidHandle<STMWaveformDigiCollection>(_conf.rawLaBrTag());
    const auto& rawdigis = *rawDigiHandle;
    _hNRAWLaBr->Fill(rawdigis.size());

    for (const auto& dg : rawdigis){
      auto const&adcs = dg.adcs();
      for (const auto& adc: adcs){
        _hRawLaBrADC->Fill(adc);
      }
      if (adcs.empty()){continue;}

      auto maxADCH = std::max_element(adcs.begin(),adcs.end());
      _hMaxRawADCLaBr -> Fill(*maxADCH);
  }
}

  if (!_conf.zsHPGeTag().empty()) {
    auto zsDigiHandle = event.getValidHandle<STMWaveformDigiCollection>(_conf.zsHPGeTag());
    const auto& zsdigis = *zsDigiHandle;
    _hNZSHPGe->Fill(zsdigis.size());

    for (const auto& dg : zsdigis){
      for (const auto& adc : dg.adcs()) _hZSHPGeADC->Fill(adc);}
  }

  if (!_conf.zsLaBrTag().empty()) {
    auto zsDigiHandle = event.getValidHandle<STMWaveformDigiCollection>(_conf.zsLaBrTag());
    const auto& zsdigis = *zsDigiHandle;
    _hNZSLaBr->Fill(zsdigis.size());

    for (const auto& dg : zsdigis){
      for (const auto& adc : dg.adcs()) _hZSLaBrADC->Fill(adc);}
  }

  if (!_conf.phHPGeTag().empty()) {
    auto phDigisHandle = event.getValidHandle<STMPHDigiCollection>(_conf.phHPGeTag());
    const auto& phDigis = *phDigisHandle;

    _hNPHHPGe->Fill(phDigis.size());

    for (const auto& phDigi: phDigis) {
      _hPulseHeightHPGe->Fill(phDigi.energy());
    }
  }

  if (!_conf.phLaBrTag().empty()) {
    auto phDigiHandle = event.getValidHandle<STMPHDigiCollection>(_conf.phLaBrTag());
    const auto& phDigis = *phDigiHandle;

    _hNPHLaBr->Fill(phDigis.size());

    for (const auto& phDigi: phDigis) {
      _hPulseHeightLaBr->Fill(phDigi.energy());
    }
  }

  if(!_conf.HPGefragSummaryTag().empty() ){
    auto summaryHandle = event.getValidHandle<STMFragmentSummaryCollection>(_conf.HPGefragSummaryTag());
    for (const auto& summary : *summaryHandle){
      _hNContainerFragsHPGe->Fill(summary.nContainerFrags());
      _hNTotalInnerFragsHPGe->Fill(summary.nInnerFrags());
      _hNGoodRawFragsHPGe ->Fill(summary.nGoodRawFrags());
      _hNGoodZSFragsHPGe ->Fill(summary.nGoodZSFrags());
      _hNGoodPHFragsHPGe ->Fill(summary.nGoodPHFrags());
      _hNZeroRawFragsHPGe ->Fill(summary.nZeroRawFrags());
      _hNZeroZSFragsHPGe ->Fill(summary.nZeroZSFrags());
      _hNZeroPHFragsHPGe ->Fill(summary.nZeroPHFrags());
      _hNEmptyRawFragsHPGe ->Fill(summary.nEmptyRawFrags());
      _hNEmptyZSFragsHPGe ->Fill(summary.nEmptyZSFrags());
      _hNEmptyPHFragsHPGe ->Fill(summary.nEmptyPHFrags());

    }
  }

  if(!_conf.LaBrfragSummaryTag().empty() ){
    auto summaryHandle = event.getValidHandle<STMFragmentSummaryCollection>(_conf.LaBrfragSummaryTag());
    for (const auto& summary : *summaryHandle){
      _hNContainerFragsLaBr->Fill(summary.nContainerFrags());
      _hNTotalInnerFragsLaBr->Fill(summary.nInnerFrags());
      _hNGoodRawFragsLaBr ->Fill(summary.nGoodRawFrags());
      _hNGoodZSFragsLaBr ->Fill(summary.nGoodZSFrags());
      _hNGoodPHFragsLaBr ->Fill(summary.nGoodPHFrags());
      _hNZeroRawFragsLaBr ->Fill(summary.nZeroRawFrags());
      _hNZeroZSFragsLaBr ->Fill(summary.nZeroZSFrags());
      _hNZeroPHFragsLaBr ->Fill(summary.nZeroPHFrags());
      _hNEmptyRawFragsLaBr ->Fill(summary.nEmptyRawFrags());
      _hNEmptyZSFragsLaBr ->Fill(summary.nEmptyZSFrags());
      _hNEmptyPHFragsLaBr ->Fill(summary.nEmptyPHFrags());
    }
  }
}

}//namespace mu2e

DEFINE_ART_MODULE(mu2e::DqmStm)
