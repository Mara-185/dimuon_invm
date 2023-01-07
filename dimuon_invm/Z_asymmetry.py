"""In this module a few functions are implemented, in order to study the
    \"Forward-Backward asymmetry" for the Z boson decay.
    The analysis is done on six different rapidity ranges of equal size and
    twelve mass bins:
    - Mass bins from article : 60 < M < 120:
        [60, 70, 78, 84, 87, 89, 91, 93, 95, 98, 104, 112, 120];
    - Rapidity bins of equal size for |y| < 2.4:
        [0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4]
"""

import ROOT
import logging
import os
import utils
import argparse

# Z Analysis

# Mass range from article : 60 < M < 120:
    # 60, 70, 78, 84, 87, 89, 91, 93, 95, 98, 104, 112, 120;
# Eta bins of equal size for |yll| < 2.4:
    # 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4

logger = logging.getLogger(__name__)

def z_analysis(infile):
    # Enable parallel analysis
    nthreads = 4
    ROOT.ROOT.EnableImplicitMT(nthreads)

    logging.info("Start the analysis...")

    # Create an RDataFrame of useful data from the root file.
    root_file = ROOT.TFile.Open(infile, "READ")
    tree_name = root_file.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_name, root_file)
    logger.debug("The RDataFrame have been created.")

    rdf2 = rdf.Define("mask", "nMuon>1").\
                Define("CutMuonPt", "Muon_pt[mask]")

    return rdf2

def mass_vs_eta(infile, rap_lim, pt_lim, bin, data_cached):
    """
    Plot the Z resonance in a mass range from 60 GeV to 120 GeV, for six
    different range of rapidity.
    """

    # RDataFrame with appropriate cuts is created.
    if data_cached == None:
        t_name = infile.replace(".root", "")
        rdf = ROOT.RDataFrame(t_name,infile)
        branch = t_name.capitalize()
        rdf_m = rdf.Filter(f"({branch}_mass>60)&&({branch}_mass<120)&&"
            f"({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]})").Filter(
            f"({branch}_pt>{pt_lim[0]})&&({branch}_pt<{pt_lim[1]})")
        del rdf
    else:
        rdf_m=data_cached.Filter(f"({branch}_mass>60)&&({branch}_mass<120)&&"
            f"({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]})").Filter(
            f"({branch}_pt>{pt_lim[0]})&&({branch}_pt<{pt_lim[1]}")
    logger.info("The RDataFrame is created.")

    # Booking histogram
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel(f"Mass in defined range rapidity",
        f"Mass in defined range rapidity;Mass [MeV];Events", bin, 60, 120),
        f"{branch}_mass")

    # Styling
    c = ROOT.TCanvas("Z mass", "Z mass")
    c.SetGrid()
    ROOT.gPad.SetLogy()
    ROOT.gStyle.SetOptStat("e")
    h.Draw()

    # Change directory to save results in "Z analysis"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}','Z analysis')))
    c.SaveAs(f"{branch}_mass_y({rap_lim[0]},{rap_lim[1]})_({pt_lim[0]},{pt_lim[1]}).png")

    # Return in main directory
    os.chdir(os.path.dirname(os. getcwd()))


def cos_vs_eta(infile, rap_lim, pt_lim, bin, data_cached):
    """
    Plot the cos of Z resonance in a mass range from 60 GeV to 120 GeV, for six
    different range of rapidity.
    """

    # RDataFrame with appropriate cuts is created.
    if data_cached == None:
        t_name = infile.replace(".root", "")
        branch = t_name.capitalize()
        rdf = ROOT.RDataFrame(t_name,infile)
        rdf_m = rdf.Filter(f"({branch}_mass>60)&&({branch}_mass<120)&&"
            f"({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]})").\
            Filter(f"({branch}_pt>{pt_lim[0]})&&({branch}_pt<{pt_lim[1]})")
        del rdf
    else:
        rdf_m=data_cached.Filter(f"({branch}_mass>60)&&({branch}_mass<120)&&"
        f"({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]})").\
        Filter(f"({branch}_pt>{pt_lim[0]})&&({branch}_pt<{pt_lim[1]})")
    logger.info("The RDataFrame is created.")

    # Booking histogram
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel(f"Cos in defined range rapidity",
        f"cos(#theta *) in defined range rapidity;Cos;Events",bin,-1,1),f"{branch}_cos")

    # Styling
    c = ROOT.TCanvas("Z cos", "Z cos")
    c.SetGrid()
    ROOT.gStyle.SetOptStat("e")
    h.Draw()

    # Change directory to save results in "Z analysis"
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}','Z analysis')))
    c.SaveAs(f"{branch}_cos_y({rap_lim[0]},{rap_lim[1]}_({pt_lim[0]},{pt_lim[1]}).png")

    # Return in main directory
    os.chdir(os.path.dirname(os. getcwd()))


def weight(infile, data_cached=None):
    """ The function computes the numerator and denominator weights for each dileptons, which
    are necessary to obtain the forward-backward asymmetry.
    """

    # RDataFrame with appropriate cuts is created.
    if data_cached == None:
        t_name = infile.replace(".root", "")
        branch = t_name.capitalize()
        rdf = ROOT.RDataFrame(t_name,infile).\
            Filter(f"({branch}_mass>60)&&({branch}_mass<120)")
    else:
        rdf=data_cached.Filter(f"({branch}_mass>60)&&({branch}_mass<120)")
    logger.info("The RDataFrame is created.")

    # Using JITting to define a function in order to calculate useful quantities
    # for further analysis.
    cppcode = """
    // The function takes in input pt, mass and cos of eah dimuon and it returns
    // a "VecOps" with useful quantities to compute the forward-backward asymmetry.
    // In particular they are "w_d" and "w_n" that are denominator and
    // numerator weights.

    //Function definition

    ROOT::VecOps::RVec<float> weights(float pt, float m, float c)
    {
        float A0, h, den, w_d, w_n;
        A0 = (pow(pt,2))/(pow(m, 2)+pow(pt,2));
        h = (0.5)*A0*(1-(3*pow(c,2)));
        den = 1+pow(c,2)+h;
        w_d = (0.5)*(pow(c,2)/pow(den,3));
        w_n = (0.5)*(abs(c)/pow(den,2));

        ROOT::VecOps::RVec<float> P{w_d,w_n};
    	return P;
    }
    """

    ROOT.gInterpreter.ProcessLine(cppcode)                                      #It compiles the C++ code
    logger.debug("The c++ function is defined.")

    rdf_w = rdf.\
        Define("w_d", f"weights({branch}_pt,{branch}_mass, {branch}_cos)[0]").\
        Define("w_n", f"weights({branch}_pt,{branch}_mass, {branch}_cos)[1]")
    logger.debug("The new columns are defined.")
    del rdf

    # A list with the names of the new columns is created in order to make a
    # snapshot of them

    # TAGLIA
    branchlist = [f"{branch}_mass", f"{branch}_pt", f"{branch}_eta", f"{branch}_phi",\
        f"{branch}_cos",f"{branch}_y","w_d", "w_n"]

    # A snapshot is done to collect the useful physiscal quantity of the dileptons
    # in a single root file.
    logger.debug("Starting snapshot...")
    rdf_w.Snapshot(f"{t_name}_w", f"{t_name}_w.root", branchlist)
    logger.info("The snapshot is done.")


def afb(infile, pt_lim):
    """ The function calculates the value of Afb (Asymmetry forward-backward)
     using the \"angular event weighting\" in the six different ranges of
     rapidity, for each mass.
     It creates different \".txt\" files with the results.
    """

    # Create RDataFrame
    t_name = infile.replace(".root", "")
    branch = t_name.replace("_w", "").capitalize()
    rdf = ROOT.RDataFrame(t_name, infile).Filter(f"({branch}_pt>{pt_lim[0]})&&"
        f"({branch}_pt<{pt_lim[1]})")
    logger.info("The RDataFrame has been created.")

    #
    for j in range(0, len(utils.RAPIDITY_BIN), 1):
        # Retrieve rapidity range bin from dictionary RAPIDITY_BIN
        y_lim = utils.RAPIDITY_BIN[f"{j}"]

        rdf_y = rdf.Filter(f"({branch}_y>{y_lim[0]})&&({branch}_y<{y_lim[1]})")
        os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Z analysis')))
        with open(f"Afb[y({y_lim[0]},{y_lim[1]})]_({pt_lim[0]}, {pt_lim[1]})"
            f"{t_name}.txt", "w+") as outf:
            print("Afb:M:A4", file=outf)
            os.chdir(os.path.dirname(os. getcwd()))
            for i in range(0, len(utils.MASS_BIN), 1):

                # Retrive mass range bin from dictonary MASS_BIN
                m_lim = utils.MASS_BIN[f"{i}"]

                rdf_m = rdf_y.Filter(f"({branch}_mass>{m_lim[0]})&&"
                    f"({branch}_mass<{m_lim[1]})")
                rdf_f = rdf_m.Filter(f"{branch}_cos>0")
                D_f = rdf_f.Sum("w_d")
                N_f = rdf_f.Sum("w_n")
                rdf_b = rdf_m.Filter(f"{branch}_cos<0")
                D_b = rdf_b.Sum("w_d")
                N_b = rdf_b.Sum("w_n")

                Afb = (3*(N_f.GetValue()-N_b.GetValue()))/(8*(D_f.GetValue()+D_b.GetValue()))
                M = m_lim[0] +((m_lim[1]-m_lim[0])/2)
                A4 = Afb*(8/3)

                print(Afb, M, A4,file=outf)
    del rdf_y, rdf_m, rdf_f, rdf_b


def afb_plot(infile):
    """ Plot the values of Afb obtain with the function \"afb\".
    """

    with open(infile) as ft:
        data = [_.rstrip('\n') for _ in ft]
    c = ROOT.TCanvas("Afb")
    file_name = infile.replace("txt", "pdf")
    c.Print(f"{file_name}[")
    for inf in data:
        t = ROOT.TTree("tree", "tree")
        t.ReadFile(inf)

        c.UseCurrentStyle()
        c.SetGrid()
        ROOT.gStyle.SetOptStat("e")
        t.Draw("Afb:M","")
        name = inf.replace("txt", "")
        c.Print(f"{file_name}", f"Title:{name}")
    c.Print(f"{file_name}]", f"Title:{name}")



if __name__ == "__main__":

    # Enable the multi-threading analysis
    nthreads = 4
    ROOT.ROOT.EnableImplicitMT(nthreads)


    # Creating the parser
    parser = argparse.ArgumentParser(description =
        "Processing the root file of data.")
    parser.add_argument("-f", "--file",required=True, type=str,
        help="The path of the nanoAOD file to analyze.")
    args = parser.parse_args()


    logger = utils.set_logger("Z analysis")

    logger.info("Starting the analysis of the root file nanoAOD...")

    rdata=z_analysis(f"{args.file}")

    rdata.Display(["Muon_pt", "CutMuonPt"])
