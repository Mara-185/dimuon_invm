"""Feynman diagram showing production of Z boson from the interaction of a quark
    and an antiquark, and its decay in two muons of opposit charge."""

import ROOT

if __name__ == "__main__":

    c = ROOT.TCanvas("Feynman diagram", "Feynman diagram")#, 10, 10, 600, 300)
    c.Range(0, 0, 120,60)
    lin = ROOT.gStyle.GetLineWidth()
    ROOT.gStyle.SetLineWidth(3)

    label = ROOT.TLatex()
    label.SetTextAlign(22)
    label.SetTextSize(.07)

    # Quarks lines
    line1 = ROOT.TLine(15, 10, 35, 30)
    line1.Draw()
    line2 = ROOT.TLine(15, 50, 35, 30)
    line2.Draw()

    label.DrawLatex(8,6,"#bar{q}")
    label.DrawLatex(8,55,"q")

    # Bosons
    gZ = ROOT.TCurlyLine(35, 30, 80, 30)
    gZ.SetWavy()
    gZ.Draw()
    label.DrawLatex(57,37.7,"#gamma/Z_{0}")

    # Muons lines
    line3 = ROOT.TLine(80, 30, 100, 10)
    line3.Draw()
    line4 = ROOT.TLine(80, 30, 100, 50)
    line4.Draw()

    label.DrawLatex(100,6,"l^{+} (#mu^{+})")
    label.DrawLatex(100,55,"l^{-} (#mu^{-})")

    label.DrawLatex(60, 5, "#bar{q}q -> #gamma/Z_{0} -> l^{+}l^{-}")

    c.Update()
    ROOT.gStyle.SetLineWidth(lin)
    c.SaveAs("Z_diagram.png")
