//
// Created by nikol on 2/21/2024.
//

#include "drawF.h"




// Method to print descriptive text to the canvas.
// (x, y) are relative coordinates (NDC).





void setUpRootStyle(){
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02, "x");
    gStyle->SetTickLength(0.015, "y");
    gStyle->SetPalette(55);
}

void drawText(double x, double y, TString txt, int align,
              double tsize) {
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

void drawParticleMarker(const Pythia8::Particle &p, int style, int col,
                        double size) {
    static auto m = new TMarker();
    m->SetMarkerStyle(style);
    m->SetMarkerSize(size);
    m->SetMarkerColor(col);
    m->DrawMarker(p.y(), p.phi());
}

//==========================================================================

// Method to draw a marker+text of a particle.

void drawParticleText(const Pythia8::Particle &p, int colourHS) {
    // Draws a marker at (y, phi) of particle. Circle for parton, star
    // for boson.
    bool isParton =  (std::abs(p.id()) <= 5 || p.id() == 21);
    int col = colourHS;
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

void drawMarker(double x, double y, int style, int col, double size) {
    auto m = new TMarker(x, y, style);
    m->SetMarkerSize(size);
    m->SetMarkerColor(col);
    m->SetNDC(true);
    m->Draw();
}

//==========================================================================