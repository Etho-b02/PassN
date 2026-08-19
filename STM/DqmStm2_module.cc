// =====================================================================
// analyzer module producing DQM histograms for Stm subsystem
// Created by Bryan Gonzalez - version 2
// =====================================================================

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
#include <algorithm>
#include <numeric>

namespace mu2e
{
    class DqmStm2 : public art::EDAnalyzer
    {
    public:
        struct Config
        {
            using Name = fhicl::Name;
            using Comment = fhicl::Comment;
            // If we are only saving PH may need to find a work around the Raw and ZS
            /*
            fhicl::Atom<art::InputTag> rawWaveformHPGeDigisTag{
                Name("rawWaveformDigisHPGeTag"), Comment("Raw STM WaveformDigi HPGe Collection"), art::InputTag()};
            fhicl::Atom<art::InputTag> rawWaveformLaBrDigisTag{
                Name("rawWaveformDigisLaBrTag"), Comment("Raw STM WaveformDigi LaBr Collection"), art::InputTag()};
            fhicl::Atom<art::InputTag> zsWaveformHPGeDigisTag{
                Name("zsWaveformDigisHPGeTag"), Comment("ZS STM WaveformDigi HPGe Collection"), art::InputTag()};
            fhicl::Atom<art::InputTag> zsWaveformLaBrDigisTag{
                Name("zsWaveformDigisLaBrTag"), Comment("ZS STM WaveformDigi LaBr Collection"), art::InputTag()};
            */

            // This Will be our main diagnostics
            fhicl::Atom<art::InputTag> phHPGeDigisTag{
                Name("phHPGeDigisTag"), Comment("STM PH Digi HPGe Collection"), art::InputTag()};
            fhicl::Atom<art::InputTag> phLaBrDigisTag{
                Name("phLaBrDigisTag"), Comment("STM PH Digi LaBr Collection"), art::InputTag()};
            fhicl::Atom<art::InputTag> fragSummaryHPGeTag{
                Name("fragSummaryHPGeTag"), Comment("STM HPGe FragmentSummary Collection"), art::InputTag()};
            fhicl::Atom<art::InputTag> fragSummaryLaBrTag{
                Name("fragSummaryLaBrTag"), Comment("STM LaBr FragmentSummary Collection"), art::InputTag()};
        };
        typedef art::EDAnalyzer::Table<Config> Parameters;

        explicit DqmStm2(const Parameters &conf);
        virtual ~DqmStm2() {}

        virtual void beginJob();
        virtual void endJob() {};
        virtual void analyze(const art::Event &e);

    private:
        Config _conf;
        // Version Information
        TH1D *_hVersion;

        // Histograms for STMWaveforms
        // Raw Waveforms
        // ZS Waveforms

        // Histograms for STM PHDigis - HPGe
        TH1D* _hNPHHPGe;          // Number of PH Digis per event
        TH1D* _hPulseHeightsHPGe; // Distribtuion of Pulse Heights, uncalibrated energy (ADC Units)
        //TH1D* _hPHTimeHPGe;       // Distribution of uncalibrated PH times, time of hit

        // Histograms for STM PHDigis - LaBr
        TH1D* _hNPHLaBr;
        TH1D* _hPulseHeightsLaBr;
        //TH1D* _hPHTimeLaBr;

        // Histograms for STM FragmentSummary --config may change
        // General Fragment information - HPGe
        //TH1D* _hNContainerFragsHPGe;
        //TH1D* _hInnerFragsHPGe;
        TH1D* _hRawFragStatusHPGe;
        TH1D* _hZSFragStatusHPGe;
        TH1D* _hPHFragStatusHPGe;

        // General Fragment information - LaBr
        //TH1D* _hNContainerFragsLaBr;
        //TH1D* _hInnerFragsLaBr;
        TH1D* _hRawFragStatusLaBr;
        TH1D* _hZSFragStatusLaBr;
        TH1D* _hPHFragStatusLaBr;
    };

    /****************************/
    DqmStm2::DqmStm2(const Parameters &conf)
        : art::EDAnalyzer(conf), _conf(conf())
    {
        // Figure out a workaround for the raw/zs or leave disabled for now
        /*
        if (!_conf.rawWaveformHPGeDigisTag().empty()){
            mayConsume<STMWaveformDigiCollection>(_conf.rawWaveformHPGeDigisTag());
        }
        if (!_conf.rawWaveformLabrDigisTag().empty()){
            mayConsume<STMWaveformDigiCollection>(_conf.rawWaveformLarBrDigisTag());
        }
        if (!_conf.zsWaveformHPGeDigisTag().empty()){
            mayConsume<STMWaveformDigiCollection>(_conf.zsWaveformHPGeDigisTag());
        }
        if (!_conf.zsWaveformLabrDigisTag().empty()){
            mayConsume<STMWaveformDigiCollection>(_conf.zsWaveformLarBrDigisTag());
        }
        */

        // Our main DQM Metrics here
        if (!_conf.phHPGeDigisTag().empty())
        {
            mayConsume<STMPHDigiCollection>(_conf.phHPGeDigisTag());
        }
        if (!_conf.phLaBrDigisTag().empty())
        {
            mayConsume<STMPHDigiCollection>(_conf.phLaBrDigisTag());
        }
        if (!_conf.fragSummaryHPGeTag().empty())
        {
            mayConsume<STMFragmentSummaryCollection>(_conf.fragSummaryHPGeTag());
        }
        if (!_conf.fragSummaryLaBrTag().empty())
        {
            mayConsume<STMFragmentSummaryCollection>(_conf.fragSummaryLaBrTag());
        }

    } // parameters

    void DqmStm2::beginJob()
    {
        art::ServiceHandle<art::TFileService> tfs;
        // Create version histogram
        _hVersion = tfs->make<TH1D>("Ver", "Version Number", 101, -0.5, 100);
        // Make histograms for our main metrics

        // Our main DQM Metrics here
        if (!_conf.phHPGeDigisTag().empty())
        {
            _hNPHHPGe = tfs->make<TH1D>("NPHDigisHPGe", "Number of PH Digis (HPGe);Number of Pulse Height Digis;Entries",
            20, -0.5,30.5);

            _hPulseHeightsHPGe = tfs->make<TH1D>("PulseHeightHPGe", "Pulse Height (HPGe);Pulse Height;Entries",
            600,-300,300);

        }

        if (!_conf.phLaBrDigisTag().empty())
        {
            _hNPHLaBr = tfs->make<TH1D>("NPHDigisLaBr", "Number of PH Digis (LaBr);Number of Pulse Height Digis;Entries",
            20, -0.5,30.5);

            _hPulseHeightsLaBr = tfs->make<TH1D>("PulseHeightLaBr", "Pulse Height (LaBr);Pulse Height;Entries",
            600,-300,300);

        }

        if (!_conf.fragSummaryHPGeTag().empty())
        {
            _hRawFragStatusHPGe = tfs->make<TH1D>(
                "RawFragStatusHPGe", "HPGe Raw Fragment Status;FragmentStatus;TotalFragments",
                4, 0.5, 4.5);
            _hRawFragStatusHPGe->GetXaxis()->SetBinLabel(1, "Prescaled");
            _hRawFragStatusHPGe->GetXaxis()->SetBinLabel(2, "Good");
            _hRawFragStatusHPGe->GetXaxis()->SetBinLabel(3, "Zero");
            _hRawFragStatusHPGe->GetXaxis()->SetBinLabel(4, "Empty");

            _hZSFragStatusHPGe = tfs->make<TH1D>(
                "ZSFragStatusHPGe", "HPGe ZS Fragment Status;FragmentStatus;TotalFragments",
                4, 0.5, 4.5);
            _hZSFragStatusHPGe->GetXaxis()->SetBinLabel(1, "Prescaled");
            _hZSFragStatusHPGe->GetXaxis()->SetBinLabel(2, "Good");
            _hZSFragStatusHPGe->GetXaxis()->SetBinLabel(3, "Zero");
            _hZSFragStatusHPGe->GetXaxis()->SetBinLabel(4, "Empty");

            _hPHFragStatusHPGe = tfs->make<TH1D>(
                "PHFragStatusHPGe", "HPGe PH Fragment Status;FragmentStatus;TotalFragments",
                3, 0.5, 3.5);
            _hPHFragStatusHPGe->GetXaxis()->SetBinLabel(1, "Good");
            _hPHFragStatusHPGe->GetXaxis()->SetBinLabel(2, "Zero");
            _hPHFragStatusHPGe->GetXaxis()->SetBinLabel(3, "Empty");
        }

        if (!_conf.fragSummaryLaBrTag().empty())
        {
            _hRawFragStatusLaBr = tfs->make<TH1D>(
                "RawFragStatusLaBr", "LaBr Raw Fragment Status;FragmentStatus;TotalFragments",
                4, 0.5, 4.5);
            _hRawFragStatusLaBr->GetXaxis()->SetBinLabel(1, "Prescaled");
            _hRawFragStatusLaBr->GetXaxis()->SetBinLabel(2, "Good");
            _hRawFragStatusLaBr->GetXaxis()->SetBinLabel(3, "Zero");
            _hRawFragStatusLaBr->GetXaxis()->SetBinLabel(4, "Empty");

            _hZSFragStatusLaBr = tfs->make<TH1D>(
                "ZSFragStatusLaBr", "LaBr ZS Fragment Status;FragmentStatus;TotalFragments",
                4, 0.5, 4.5);
            _hZSFragStatusLaBr->GetXaxis()->SetBinLabel(1, "Prescaled");
            _hZSFragStatusLaBr->GetXaxis()->SetBinLabel(2, "Good");
            _hZSFragStatusLaBr->GetXaxis()->SetBinLabel(3, "Zero");
            _hZSFragStatusLaBr->GetXaxis()->SetBinLabel(4, "Empty");

            _hPHFragStatusLaBr = tfs->make<TH1D>(
                "PHFragStatusLaBr", "LaBr PH Fragment Status;FragmentStatus;TotalFragments",
                3, 0.5, 3.5);
            _hPHFragStatusLaBr->GetXaxis()->SetBinLabel(1, "Good");
            _hPHFragStatusLaBr->GetXaxis()->SetBinLabel(2, "Zero");
            _hPHFragStatusLaBr->GetXaxis()->SetBinLabel(3, "Empty");
        }

    } // beginJob

    /******************************/
    void DqmStm2::analyze(const art::Event &event)
    {
        _hVersion->Fill(0.0);

        // Our main DQM Metrics here

        // Fill ph HPGe Digis here
        if (!_conf.phHPGeDigisTag().empty())
        {
            auto phDigiHandle = event.getValidHandle<STMPHDigiCollection>(_conf.phHPGeDigisTag());
            const auto& phDigis = *phDigiHandle;

            _hNPHHPGe->Fill(phDigis.size());

            for (const auto& phDigi : phDigis) {
                _hPulseHeightsHPGe->Fill(phDigi.energy());
            }

        }

        // Fill ph LaBR Digis here
        if (!_conf.phLaBrDigisTag().empty())
        {
            auto phDigiHandle = event.getValidHandle<STMPHDigiCollection>(_conf.phLaBrDigisTag());
            const auto& phDigis = *phDigiHandle;

            _hNPHLaBr->Fill(phDigis.size());

            for (const auto& phDigi : phDigis) {
                _hPulseHeightsLaBr->Fill(phDigi.energy());
            }
        }

        // Fill HPGe Frag status
        if (!_conf.fragSummaryHPGeTag().empty())
        {
            auto summaryHandle = event.getValidHandle<STMFragmentSummaryCollection>(_conf.fragSummaryHPGeTag());
            for (const auto &summary : *summaryHandle)
            {
                _hRawFragStatusHPGe->Fill(1, summary.nPrescaledRawFrags());
                _hRawFragStatusHPGe->Fill(2, summary.nGoodRawFrags());
                _hRawFragStatusHPGe->Fill(3, summary.nZeroRawFrags());
                _hRawFragStatusHPGe->Fill(4, summary.nEmptyRawFrags());

                _hZSFragStatusHPGe->Fill(1, summary.nPrescaledZSFrags());
                _hZSFragStatusHPGe->Fill(2, summary.nGoodZSFrags());
                _hZSFragStatusHPGe->Fill(3, summary.nZeroZSFrags());
                _hZSFragStatusHPGe->Fill(4, summary.nEmptyZSFrags());

                _hPHFragStatusHPGe->Fill(1, summary.nGoodPHFrags());
                _hPHFragStatusHPGe->Fill(2, summary.nZeroPHFrags());
                _hPHFragStatusHPGe->Fill(3, summary.nEmptyPHFrags());
            }
        }

        // Fill LaBr Frag Status
        if (!_conf.fragSummaryLaBrTag().empty())
        {
            auto summaryHandle = event.getValidHandle<STMFragmentSummaryCollection>(_conf.fragSummaryLaBrTag());
            for (const auto &summary : *summaryHandle)
            {
                _hRawFragStatusLaBr->Fill(1, summary.nPrescaledRawFrags());
                _hRawFragStatusLaBr->Fill(2, summary.nGoodRawFrags());
                _hRawFragStatusLaBr->Fill(3, summary.nZeroRawFrags());
                _hRawFragStatusLaBr->Fill(4, summary.nEmptyRawFrags());

                _hZSFragStatusLaBr->Fill(1, summary.nPrescaledZSFrags());
                _hZSFragStatusLaBr->Fill(2, summary.nGoodZSFrags());
                _hZSFragStatusLaBr->Fill(3, summary.nZeroZSFrags());
                _hZSFragStatusLaBr->Fill(4, summary.nEmptyZSFrags());

                _hPHFragStatusLaBr->Fill(1, summary.nGoodPHFrags());
                _hPHFragStatusLaBr->Fill(2, summary.nZeroPHFrags());
                _hPHFragStatusLaBr->Fill(3, summary.nEmptyPHFrags());
            }
        }

    } // analyze

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::DqmStm2)
