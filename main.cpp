#include <iostream>
#include <string>

#include "Pythia8/Pythia.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH2D.h"
#include "TMath.h"
#include "TPave.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "drawF.h"

int main() {


    Pythia8::Pythia pythia;
    pythia.readFile("../config1.cmnd");
    pythia.init();

    setUpRootStyle();
    auto canvas = new TCanvas();
    canvas->SetMargin(0.06, 0.02, 0.08, 0.06);
    auto pTflow = createTH2D();




    int colHS = kBlack, colPos = kRed, colNeg = kBlue;
    int colNeut = kGreen + 3, colPU = kGray + 1;
    TString pdf = "result.pdf";
    canvas->Print(pdf + "[");

    std::map<TString, fastjet::JetDefinition> jetDefs;


    double R = 0.3;
    bool doK_t = true, doAntiK_t = true, doCambridgeAachen = false;
    int pTmin_jet = 25;
    double pTmin_hadron = 1, yMax = 4;


    if (doAntiK_t)
        jetDefs["Anti-#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
                fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    if (doK_t)
        jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
                fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    if (doCambridgeAachen)
        jetDefs["Cambridge-Aachen jets,  #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
                fastjet::cambridge_algorithm, R, fastjet::E_scheme, fastjet::Best);

    auto &event = pythia.event;

    for (int iEvent = 0; iEvent < pythia.mode("Main:numberOfEvents"); ++iEvent) { //choosing final particles only
        pTflow->Reset();
        if (!pythia.next()) continue;
        std::vector<Pythia8::Particle> particles_histogram;
        std::vector<fastjet::PseudoJet> stable_particles;
        for (int i = 0; i < event.size(); ++i) {
            auto &p = event[i];
            if (not p.isFinal()) continue;
            stable_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
            particles_histogram.push_back(p);
        }

        //Ghost are needed otherwise jet images suck
        fastjet::PseudoJet ghost;
        double pTghost = 1e-100;
        for (int iy= 1;iy<= 400/2; ++iy) {
            for (int iphi= 1;iphi<= 314/2; ++iphi) {
                double y   = pTflow->GetXaxis()->GetBinCenter(iy);
                double phi = pTflow->GetYaxis()->GetBinCenter(iphi);
                ghost.reset_momentum_PtYPhiM(pTghost, y, phi, 0);
                stable_particles.push_back(ghost);
            }
        }


        canvas->SetLogz(); //log the z axis, so jets are more clearly seen
        canvas->SetRightMargin(0.14);



        for (auto jetDef: jetDefs) {
            pTflow->Reset();


            fastjet::ClusterSequence clustSeq(stable_particles, jetDef.second);
            auto jets = sorted_by_pt(clustSeq.inclusive_jets(pTmin_jet));
            // Fill the pT flow.
            // For each jet:

            for (auto jet: jets) {
                // For each particle:
                for (auto c: jet.constituents()) {
                    if (c.pt() > 1e-50) continue; //gets rid of bubbles in jets
                    pTflow->Fill(c.rap(), c.phi_std(), jet.pt());
                }
            }

            pTflow->GetZaxis()->SetRangeUser(pTmin_jet / 4, pTflow->GetBinContent(pTflow->GetMaximumBin()) * 4);
            pTflow->GetZaxis()->SetMoreLogLabels();
            pTflow->Draw("colz");

            drawParticles_histogram(particles_histogram);

            drawdrawLegend();
            canvas->Print(pdf);
        }
    }
    canvas->Print(pdf + "]");
    printf("Produced %s\n\n", pdf.Data());

    delete pTflow;
    delete canvas;

    return 0;
}
