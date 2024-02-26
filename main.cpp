#include <iostream>
#include <string>
#include "Pythia8/Pythia.h"

#include "drawF.h"

int main() {
    Pythia8::Pythia pythia;
    pythia.readFile("../config1.cmnd");
    pythia.init();

    setUpRootStyle();
    auto can = new TCanvas();
    can->SetMargin(0.06, 0.02, 0.08, 0.06);
    auto pTflow = new TH2D("", ";Rapidity #it{y};Azimuth #it{#phi};Jet #it{p}_{T} [GeV]",
                           400 / 2, -4, 4, 314 / 2, -TMath::Pi(), TMath::Pi());

    pTflow->GetYaxis()->SetTitleOffset(0.8);
    pTflow->GetZaxis()->SetTitleOffset(1.1);

    std::map<TString, fastjet::JetDefinition> jetDefs;


    double R = 0.3;
    bool doK_t = true, doAntiK_t = true, doCambridgeAachen = false;


    if (doAntiK_t)
        jetDefs["Anti-#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
                fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    if (doAntiK_t)
        jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
                fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    if (doCambridgeAachen)
        jetDefs["Cambridge-Aachen jets,  #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
                fastjet::cambridge_algorithm, R, fastjet::E_scheme, fastjet::Best);

    auto &event = pythia.event;

    for (int iEvent = 0; iEvent < pythia.mode("Main:numberOfEvents"); ++iEvent) { //choosing final particles only
        if (!pythia.next()) continue;
        std::vector<Pythia8::Particle> particles_histogram;
        std::vector<fastjet::PseudoJet> stable_particles;
        for (int i = 0; i < event.size(); ++i) {
            auto &p = event[i];
            if (not p.isFinal()) continue;
            stable_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
            particles_histogram.push_back(p);
        }

        for (auto jetDef: jetDefs) {
            fastjet::ClusterSequence clustSeq(stable_particles, jetDef.second);
            auto jets = sorted_by_pt(clustSeq.inclusive_jets());
            // Fill the pT flow.
            pTflow->Reset();
            // For each jet:
            for (auto jet: jets) {
                // For each particle:
                for (auto c: jet.constituents()) {
                    if (c.pt() > 1e-50) continue;
                    pTflow->Fill(c.rap(), c.phi_std(), jet.pt());
                }
            }
        }

        for (auto &p: particles_histogram) {

            if (p.charge() > 0) {
                drawParticleMarker(p, 5, colPos, 0.8);
            } else if (p.charge() < 0) {
                drawParticleMarker(p, 5, colNeg, 0.8);
            } else {
                drawParticleMarker(p, 21, colNeut, 0.4);
                drawParticleMarker(p, 5, colNeut, 0.8);
            }

        }


    }
    return 0;
}
