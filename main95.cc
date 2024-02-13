// main95.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Dag Gillberg <dag.gillberg@cern.ch>

// Keywords: root; jets; event display

// This is a program to use ROOT to visualize different jet algoritms.
// The produced figure is used in the article "50 years of Quantum
// Chromodynamics" in celebration of the 50th anniversary of QCD (EPJC).

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

//==========================================================================

// Hard-coded settings

// Jet and hadron pT thresholds.
// Will only show particles with pT > pTmin and |y| < yMax.
double pTmin_jet = 25;
double pTmin_hadron = 1;
double yMax = 4;

// Amount of pileup. Average number of inelastic pp collisions per event
// (=bunch-crossing). Set to zero to turn off pileup.
double mu = 60;

// Style format. Colours used by various drawn markers.
int colHS = kBlack, colPos = kRed, colNeg = kBlue;
int colNeut = kGreen + 3, colPU = kGray + 1;

using namespace Pythia8;

//==========================================================================

// Method to print descriptive text to the canvas.
// (x, y) are relative coordinates (NDC).

void drawText(double x, double y, TString txt, int align= 11,
  double tsize= 0.032) {
  static auto tex = new TLatex();
  tex->SetTextAlign(align);
  tex->SetTextSize(tsize);
  tex->SetTextFont(42);
  tex->SetNDC();
  tex->DrawLatex(x, y, txt);
}

//==========================================================================

// Text to draw a marker at the (y, phi) coordinates of a particle.
// Absolute coordinates.

void drawParticleMarker(const Particle &p, int style, int col,
  double size= 1.0) {
  static auto m = new TMarker();
  m->SetMarkerStyle(style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(col);
  m->DrawMarker(p.y(), p.phi());
}

//==========================================================================

// Method to draw a marker+text of a particle.

void drawParticleText(const Particle &p) {
  // Draws a marker at (y, phi) of particle. Circle for parton, star
  // for boson.
  bool isParton =  (std::abs(p.id()) <= 5 || p.id() == 21);
  int col = colHS;
  drawParticleMarker( p, isParton?20:29, col, isParton?0.8:1.2);

  // Format the name-string of the particle according to ROOT's TLatex.
  // Print the text right under the marker.
  TString name = p.name();
  if (name.Contains("bar")) name = "#bar{" + name.ReplaceAll("bar", "") + "}";
  name.ReplaceAll("+", "^{+}").ReplaceAll("-", "^{-}").ReplaceAll("h0", "H");
  static auto tex = new TLatex();
  tex->SetTextSize(0.03);
  tex->SetTextFont(42);
  tex->SetTextAlign(11);
  tex->SetTextColor(col);
  tex->DrawLatex(p.y() + 0.1, p.phi() - 0.1, "#it{" + name + "}");
}

//==========================================================================

// Draws a box for text to appear.

void drawLegendBox(double x1, double y1, double x2, double y2) {
  static auto *box = new TPave(x1, y1, x2, y2, 1, "ndc");
  box->SetFillColor(kWhite);
  box->Draw();
}

//==========================================================================

// Draw a marker for legend.

void drawMarker(double x, double y, int style, int col, double size= 1.0) {
  auto m = new TMarker(x, y, style);
  m->SetMarkerSize(size);
  m->SetMarkerColor(col);
  m->SetNDC(true);
  m->Draw();
}

//==========================================================================

// Example main program to vizualize jet algorithms.

int main() {
	Pythia pythia;
	Event& event = pythia.event;
	pythia.readFile("main95conf.cmnd");
	
	
	
	
	
	
	
	// Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

  // Book histograms.
  Hist pThard("process pT scale", 100, 0., 500.);
  Hist nFinal("final particle multiplicity", 100, -0.5, 1599.5);
  Hist nCharged("charged particle multiplicity", 100, -0.5, 799.5);
  Hist dndy("dn/dy for charged particles", 100, -10., 10.);
  Hist dndeta("dn/d(eta) for charged particles", 100, -10., 10.);
  Hist dndpT("dn/dpT for charged particles", 100, 0., 10.);
  Hist epCons("deviation from energy-momentum conservation", 100, 0., 1e-6);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Fill hard scale of event.
    pThard.fill( pythia.info.pTHat() );

    // Loop over final particles in the event.
    int  nFin = 0;
    int  nChg = 0;
    Vec4 pSum;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

       // Analyze all particles.
      nFin++;
      pSum += event[i].p();

      // Analyze charged particles and fill histograms.
      if (event[i].isCharged()) {
        ++nChg;
        dndy.fill( event[i].y() );
        dndeta.fill( event[i].eta() );
        dndpT.fill( event[i].pT() );
      }

    // End of particle loop. Fill global properties.
    }
    nFinal.fill( nFin );
    nCharged.fill( nChg );
    pSum /= event[0].e();
    double epDev = abs(pSum.e() - 1.) + abs(pSum.px()) + abs(pSum.py())
      + abs(pSum.pz());
    epCons.fill(epDev);

  // End of event loop.
  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();
  double sigma = pythia.info.sigmaGen();
  pThard   *= sigma * 1e6 * 0.2 / nEvent;
  nFinal   *= 1. / (16. * nEvent);
  nCharged *= 1. / (8. * nEvent);
  dndy     *=  5. / nEvent;
  dndeta   *=  5. / nEvent;
  dndpT    *= 10. / nEvent;
  cout << pThard << nFinal << nCharged << dndy << dndeta << dndpT << epCons;

  // Write Python code that can generate a PDF file with the spectra.
  // Assuming you have Python installed on your platform, do as follows.
  // After the program has run, type "python main03plot.py" (without the " ")
  // in a terminal window, and open "out03plot.pdf" in a PDF viewer.
  // Colours and other choices can be omitted, but are shown as illustration.
  HistPlot hpl("main03plot");
  hpl.plotFrame("out03plot", pThard, "$p_{\\perp}$ scale of hard interaction",
    "$p_{\\perp}$ (GeV)",
    "$\\mathrm{d}\\sigma/\\mathrm{d}p_{\\perp}$ (nb/GeV)",
    "h", "$p_{\\perp}$ of $2 \\to 2$ process", true);
  hpl.frame("", "Total and charged particle multiplicities",
    "$n$", "$\\mathrm{d}P/\\mathrm{d}n$");
  hpl.add( nFinal, "h,royalblue", "total");
  hpl.add( nCharged, "h,orange", "charged (even only!)");
  hpl.plot();
  hpl.frame( "", "Charged (pseudo)rapidity distribution", "$y$ or $\\eta$",
    "$\\mathrm{d}n_{\\mathrm{charged}}/\\mathrm{d}(y/\\eta)$");
  hpl.add( dndy, "-", "$\\mathrm{d}n_{\\mathrm{charged}}/\\mathrm{d}y$");
  hpl.add( dndeta, "--,magenta",
    "$\\mathrm{d}n_{\\mathrm{charged}}/\\mathrm{d}\\eta$");
  hpl.plot();
  hpl.plotFrame("", dndpT, "Charged $p_{\\perp}$ spectrum",
    "$p_{\\perp}$ (GeV)", "$\\mathrm{d}n_{\\mathrm{charged}}/\\mathrm{d}"
    "p_{\\perp}$ (GeV$^{-1}$)", "", "charged", true);

  // Done.

	return 0;
}
