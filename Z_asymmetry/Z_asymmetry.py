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
import time
from pathlib import Path


# Z Analysis

# Mass range from article : 60 < M < 120:
    # 60, 70, 78, 84, 87, 89, 91, 93, 95, 98, 104, 112, 120;
# Eta bins of equal size for |yll| < 2.4:
    # 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4

logger = logging.getLogger(__name__)

def retrieve_dataset(findex, start, stop):
    """
    The function takes as argument the file index with all name of root files
    of the chosen dataset and the two extremities of the file's range to analyze.

    :param findex : index of root files which make up the entire dataset
    :type findex: txt file, required.
    :param start: first index of the range of files to analyze
    :type start: int, required
    :param stop: second index of the range of files to analyze
    :type stop: int, required
    """

    with open(findex) as inf:
        all_data = [_.rstrip('\n') for _ in inf]
        logger.info(f"Total number of data files: {len(all_data)}\n")
        data = all_data[start:stop]
        times = []

        root_files = ROOT.std.vector("string")()
        for i,d in enumerate(data):
            logger.info(f"Analyzing file number {i+1} :\n{d}")
            start1 = time.time()
            root_files.push_back(z_main(d, i))
            stop1 = time.time()
            times.append(stop1-start1)
            logger.info(f"Time : {stop1-start1}")

    return root_files, times


def z_analysis(infile, iteration):
    """
    The function create an RDataFrame to do the right cuts on data and it creates
    a snapshot of the events that passed the cuts.
    ????????????????????????????????????????????????????????

        - infile: a vector of strings which contains the files' names to analyze.
            This is the output of the function \"retrieve_dataset\".
        - iteration: the index of the file which is been analyzing ???

    :param infile: vector of strings with the name of root files
    :type infile: vector, required
    :param iteration: number of file analyzed
    :type iteration: int, required
    """

    logging.info("Start the analysis...")

    # Create an RDataFrame of useful data from the root file.
    root_file = ROOT.TFile.Open(infile, "READ")
    tree_name = root_file.GetListOfKeys().At(0).GetName()
    rdf = ROOT.RDataFrame(tree_name, root_file)
    logger.debug("The RDataFrame have been created.")

    rdf_cut=rdf.Filter("nMuon>1","Selection on events with at least two muons").\
        Define("mask", "(Muon_isTracker==1)&&(Muon_isGlobal==1)&&"
        "(Muon_gChi2<10)&&(Muon_dxyBest<0.2)&&(Muon_dxyBest>-0.2)&&"
        "(Muon_pfRelIso03_all<0.1)&&(Muon_eta>-2.4)&&(Muon_eta<2.4)").\
        Define("Mu_pt", "Muon_pt[mask]").\
        Define("Mu_pterr", "Muon_ptErr[mask]").\
        Define("Mu_mass", "Muon_mass[mask]").\
        Define("Mu_charge", "Muon_charge[mask]").\
        Define("Mu_phi", "Muon_phi[mask]").\
        Define("Mu_eta", "Muon_eta[mask]").\
        Filter("(Mu_pt.size())>=nMuon").\
        Filter("(Mu_pt[0]>25)&&(Mu_pt[1])>15",
               "Trigger on the leading and the next leading muons").\
        Filter("Mu_charge[0]!=Mu_charge[1]",
               "Selection on muons with opposite charge")

    if iteration==0:

        ROOT.gInterpreter.Declare("""
        ROOT::VecOps::RVec<double> dilepton_vec(float pt0, float eta0, float phi0,
            float mass0, float pt1, float eta1, float phi1, float mass1)
        {
        	ROOT::Math::PtEtaPhiMVector p0, p1;
        	p0.SetCoordinates(pt0, eta0, phi0, mass0);
        	p1.SetCoordinates(pt1, eta1, phi1, mass1);
        	ROOT::VecOps::RVec<double> P{(p1+p0).Pt(), (p1+p0).Eta(), (p1+p0).Phi(),
                    (p1+p0).M()};
        	return P;
        }

        ROOT::VecOps::RVec<float> cos_rapidity(float pt0, float eta0, float phi0,
                float mass0, float pt1, float eta1, float phi1, float mass1)
        {
            float pz0, pz1, E0, E1, P0_1, P0_2, P1_1, P1_2, cos, numer, denom, mll,
                ptt, pzll, y;
            ROOT::Math::PtEtaPhiMVector p0, p1;
      		p0.SetCoordinates(pt0, eta0, phi0, mass0);
      		p1.SetCoordinates(pt1, eta1, phi1, mass1);

            pz0 = pt0*sinh(eta0);
            E0 = sqrt(pow(pz0,2)+pow(pt0,2)+pow(mass0,2));
            pz1 = pt1*sinh(eta1);
            E1 = sqrt(pow(pz1,2)+pow(pt1,2)+pow(mass1,2));
            P0_1 = (E0+pz0)/sqrt(2);
            P0_2 = (E0-pz0)/sqrt(2);
            P1_1 = (E1+pz1)/sqrt(2);
            P1_2 = (E1-pz1)/sqrt(2);
            numer = 2*((P0_1*P1_2) - (P0_2*P1_1));
            mll = pow((p1+p0).M(),2);
            ptt = pow((p1+p0).Pt(),2);
            denom = sqrt(mll*(mll+ptt));
            pzll = pz0+pz1;

            cos = (numer/denom)*(pzll/abs(pzll));
            y = (0.5)*log(((E0+E1)+(pz0+pz1))/((E0+E1)-(pz0+pz1)));

            ROOT::VecOps::RVec<float> A{cos,y};
      		return A;
      	}
        """)

        logger.debug("The two c++ functions have been defined.")

    rdf_dimu = rdf_cut.\
        Define("Dimuon_mass", \
            "dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[3]").\
        Define("Dimuon_pt",\
            "dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[0]").\
        Define("Dimuon_eta", \
            "dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[1]").\
        Define("Dimuon_phi",\
            "dilepton_vec(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[2]").\
        Define("Dimuon_cos",\
            "cos_rapidity(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[0]").\
        Define("Dimuon_y",\
            "cos_rapidity(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_mass[0], \
             Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_mass[1])[1]")

    del rdf

    ROOT.RDF.SaveGraph(rdf_dimu, "rdf_mu.dot")
    rdf_dimu.Report().Print()

    branchlist_mu = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi",
        "Dimuon_cos", "Dimuon_y"]
    rdf_dimu.Snapshot("dimuon", "dimuon.root", branchlist_mu)
    logger.info("The snapshot is done.")
    root_file.Close()

def dimu(infile, n):
    pass
    #####################################################
    ######ARE THEY SORT BY THE GREATEST IN PT?###########


    # f=ROOT.TFile.Open(infile, "READ")
    # tree_name = f.GetListOfKeys().At(0).GetName()
    # t = ROOT.TTree("tree", "tree")
    # t = f.Get(tree_name)

    # t.Scan("Mu_pt:Mu_pterr:Mu_mass:Mu_charge:Mu_phi:Mu_eta:nMuon")

def mass_vs_eta(rdf0, rap_lim, pt_lim, bin, data_cached=None, infile=None):
    """
    Plot the Z resonance in a mass range from 60 GeV to 120 GeV, for six
    different range of rapidity.
    """
    branch="Dimuon"

    # RDataFrame with appropriate cuts is created.
    if data_cached == None:
        if infile == None:
            rdf_m=rdf0.Filter(f"({branch}_mass>60)&&({branch}_mass<120)").Filter(
                f"(({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]}))||(({branch}_y>-{rap_lim[1]})&&({branch}_y<-{rap_lim[0]}))", "Cut on mass and rapidity").Filter(
                f"({branch}_pt>{pt_lim[0]})&&({branch}_pt<{pt_lim[1]})", "Cut on pt")
            rdf_m.Report().Print()
        else:
            t_name = infile.replace(".root", "")
            rdf = ROOT.RDataFrame(t_name,infile)
            branch = t_name.capitalize()
            rdf_m = rdf.Filter(f"({branch}_mass>60)&&({branch}_mass<120)").Filter(
                f"(({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]}))||(({branch}_y>-{rap_lim[1]})&&({branch}_y<-{rap_lim[0]}))").Filter(
                f"({branch}_pt>{pt_lim[0]})&&({branch}_pt<{pt_lim[1]})")
            del rdf
    else:
        rdf_m=data_cached.Filter(f"({branch}_mass>60)&&({branch}_mass<120)").Filter(
            f"(({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]}))||(({branch}_y>-{rap_lim[1]})&&({branch}_y<-{rap_lim[0]}))").Filter(
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

def cos_vs_eta(rdf0, rap_lim, pt_lim, bin, data_cached=None, infile=None):
    """
    Plot the cos of Z resonance in a mass range from 60 GeV to 120 GeV, for six
    different range of rapidity.
    """

    branch="Dimuon"

    # RDataFrame with appropriate cuts is created.
    if data_cached == None:
        if infile == None:
            rdf_m = rdf0.Filter(f"({branch}_mass>60)&&({branch}_mass<120)").Filter(
                f"(({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]}))||(({branch}_y>-{rap_lim[1]})&&({branch}_y<-{rap_lim[0]}))", "Cut on mass and rapidity").\
                Filter(f"({branch}_pt>{pt_lim[0]})&&({branch}_pt<{pt_lim[1]})", "Cut on pt")
            rdf_m.Report().Print()
        else:
            t_name = infile.replace(".root", "")
            branch = t_name.capitalize()
            rdf = ROOT.RDataFrame(t_name,infile)
            rdf_m = rdf.Filter(f"({branch}_mass>60)&&({branch}_mass<120)").Filter(
                f"(({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]}))||(({branch}_y>-{rap_lim[1]})&&({branch}_y<-{rap_lim[0]}))").\
                Filter(f"({branch}_pt>{pt_lim[0]})&&({branch}_pt<{pt_lim[1]})")
            del rdf
    else:
        rdf_m=data_cached.Filter(f"({branch}_mass>60)&&({branch}_mass<120)").Filter(
        f"(({branch}_y>{rap_lim[0]})&&({branch}_y<{rap_lim[1]}))||(({branch}_y>-{rap_lim[1]})&&({branch}_y<-{rap_lim[0]}))").\
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

def weight(infile, n, data_cached=None):
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
        rdf=data_cached.Filter(f"({branch}_mass>60)&&({branch}_mass<120)", "Cut on Z resoanance")
    logger.info("The RDataFrame is created.")

    rdf.Report().Print()

    if n==0:
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
    branchlist = [f"{branch}_mass", f"{branch}_pt", f"{branch}_eta", f"{branch}_phi",\
        f"{branch}_cos",f"{branch}_y","w_d", "w_n"]

    # A snapshot is done to collect the useful physiscal quantity of the dileptons
    # in a single root file.
    logger.debug("Starting snapshot...")
    rdf_w.Snapshot(f"{t_name}_w", f"{t_name}_w{n}.root", branchlist)
    logger.info("The snapshot is done.")

    return str(f"{t_name}_w{n}.root")


def afb(rdf0, pt_lim, infile=None):
    """ The function calculates the value of Afb (Asymmetry forward-backward)
     using the \"angular event weighting\" in the six different ranges of
     rapidity, for each mass.
     It creates different \".txt\" files with the results.
    """
    branch="Dimuon"
    t_name="dimuon_w"

    # Create RDataFrame
    if infile==None:
        rdf=rdf0.Filter(f"({branch}_pt>{pt_lim[0]})&&"
            f"({branch}_pt<{pt_lim[1]})")
        rdf.Report().Print()
        logger.info("The RDataFrame has been created.")
    else:
        t_name = infile.replace(".root", "")
        branch = t_name.replace("_w", "").capitalize()
        rdf = ROOT.RDataFrame(t_name, infile).Filter(f"({branch}_pt>{pt_lim[0]})&&"
            f"({branch}_pt<{pt_lim[1]})")
        logger.info("The RDataFrame has been created.")

    #
    for j in range(0, len(utils.RAPIDITY_BIN), 1):
        # Retrieve rapidity range bin from dictionary RAPIDITY_BIN
        y_lim = utils.RAPIDITY_BIN[f"{j}"]

        rdf_y = rdf.Filter(f"(({branch}_y>{y_lim[0]})&&({branch}_y<{y_lim[1]}))||(({branch}_y>-{y_lim[1]})&&({branch}_y<-{y_lim[0]}))")
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

                #######################GESTISCI ZeroDivisionError
                try:
                    Afb = (3*(N_f.GetValue()-N_b.GetValue()))/(8*(D_f.GetValue()+D_b.GetValue()))
                    A4 = Afb*(8/3)
                except ZeroDivisionError:
                    Afb = -100
                    print(Afb)
                    A4 = -100

                M = m_lim[0] +((m_lim[1]-m_lim[0])/2)

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
        t.Draw("Afb:M","Afb!=-100")
        name = inf.replace("txt", "")
        c.Print(f"{file_name}", f"Title:{name}")
    c.Print(f"{file_name}]")#, f"Title:{name}")

def z_main(url, n):
    """
    """
    # Selection on events
    z_analysis(url, n)

    # AFB (Forward_Backward Asymmetry)
    data_string = weight("dimuon.root", n)
    return data_string


if __name__ == "__main__":

    # Set ROOT environment
    ROOT.gROOT.SetBatch()
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.RooMsgService.instance().setSilentMode(True)

    # Style of the canvas
    myStyle = ROOT.TStyle('My Style', 'My Style')
    myStyle.SetCanvasColor(0)
    myStyle.SetMarkerSize(1.2)
    myStyle.SetMarkerColor(861)
    myStyle.SetMarkerStyle(20)
    myStyle.SetHistLineColor(9)
    myStyle.SetLabelFont(42)
    myStyle.SetTitleFont(22)
    myStyle.SetLabelSize(0.035)
    ROOT.gROOT.SetStyle("My Style")
    myStyle.cd()

    # Start the the timer
    timer = ROOT.TStopwatch()
    timer.Start()

    start =time.time()

    # Creating the parser
    parser = argparse.ArgumentParser(description =
        "Processing the root file of data.")
    parser.add_argument("-f", "--file",required=True, type=str,
        help="The name of the file index to pass to the script. "
             "It has to contain the path of the different root files to analyze,"
             "which are part of the dataset.E.g.:\"Run2011A_SingleMu_index.txt\"")
    parser.add_argument("-d", "--display", required=True, type=bool,
        help="Show RDataFrame")
    parser.add_argument("-l", "--limit", nargs="+",type=int,
        help="Range of files to analyze, taken from the file index."
             " Format: -l start stop. \"stop\" is excuded.")
    parser.add_argument("-a", "--analysis", required=True, type=bool,
        help="If False, retrieve data from internet, otherwise not.")
    args = parser.parse_args()

    ######## CONTROLLA INPUT

    # Enable the multi-threading analysis
    if not args.display:
        nthreads =10
        ROOT.ROOT.EnableImplicitMT(nthreads)

    logger = utils.set_logger("Z analysis")

    # Z ANALYSIS
    logger.info("Starting the analysis of the Forward-Backward asymmetry"
                " of Z resonance.")

    # A new folder to collect plots and data is created
    os.makedirs("Z analysis", exist_ok=True)
    logger.debug("The new directory \"Z analysis\" is created")

    try:
        limits = tuple(args.limit)
    except ValueError:
        logging.error("The inserted values are not right. Check on signature "
                      "in the help's parser.")

    if args.analysis==False:
        root_files, times = retrieve_dataset(args.file, limits[0], limits[1])
    else:
        root_files = ROOT.std.vector("string")()
        root_files2 = ["dimuon_w0.root", "dimuon_w1.root", "dimuon_w2.root", "dimuon_w3.root", "dimuon_w4.root",
         "dimuon_w5.root", "dimuon_w6.root", "dimuon_w7.root", "dimuon_w8.root", "dimuon_w9.root", "dimuon_w10.root",
          "dimuon_w11.root", "dimuon_w12.root", "dimuon_w13.root", "dimuon_w14.root", "dimuon_w15.root", "dimuon_w16.root",
          "dimuon_w17.root", "dimuon_w18.root", "dimuon_w19.root", "dimuon_w20.root", "dimuon_w21.root", "dimuon_w22.root",
          "dimuon_w23.root", "dimuon_w24.root", "dimuon_w25.root", "dimuon_w26.root", "dimuon_w27.root", "dimuon_w28.root",
          "dimuon_w29.root", "dimuon_w30.root", "dimuon_w31.root", "dimuon_w32.root", "dimuon_w33.root", "dimuon_w34.root",
          "dimuon_w35.root", "dimuon_w36.root", "dimuon_w37.root", "dimuon_w38.root", "dimuon_w39.root", "dimuon_w40.root",
          "dimuon_w41.root", "dimuon_w42.root", "dimuon_w43.root", "dimuon_w44.root", "dimuon_w45.root", "dimuon_w46.root",
          "dimuon_w47.root", "dimuon_w48.root", "dimuon_w49.root", "dimuon_w50.root", "dimuon_w51.root", "dimuon_w52.root",
          "dimuon_w53.root", "dimuon_w54.root", "dimuon_w55.root", "dimuon_w56.root", "dimuon_w57.root", "dimuon_w58.root",
          "dimuon_w59.root", "dimuon_w60.root", "dimuon_w61.root", "dimuon_w62.root", "dimuon_w63.root", "dimuon_w64.root",
          "dimuon_w65.root", "dimuon_w66.root", "dimuon_w67.root", "dimuon_w68.root", "dimuon_w69.root", "dimuon_w70.root",
          "dimuon_w71.root", "dimuon_w72.root", "dimuon_w73.root", "dimuon_w74.root", "dimuon_w75.root", "dimuon_w76.root",
          "dimuon_w77.root", "dimuon_w78.root", "dimuon_w79.root", "dimuon_w80.root", "dimuon_w81.root", "dimuon_w82.root",
          "dimuon_w83.root", "dimuon_w84.root", "dimuon_w85.root", "dimuon_w86.root", "dimuon_w87.root", "dimuon_w88.root",
          "dimuon_w89.root", "dimuon_w90.root", "dimuon_w91.root", "dimuon_w92.root", "dimuon_w93.root", "dimuon_w94.root",
          "dimuon_w95.root", "dimuon_w96.root", "dimuon_w97.root", "dimuon_w98.root"]

        for t in root_files2:
            root_files.push_back(t)

    # with open(f"times.txt", "w+") as outf:
    #     for t in times:
    #         print(t, file=outf)

    all_data = ROOT.RDataFrame("dimuon_w", root_files)
    N = all_data.Count().GetValue()
    logger.info(f"Total number of events that passed the cuts: {N}.")


    # Plot the distributions of mass and cos(theta*) in different range of
    # rapidity.
    pt_lim = (10,100)
    #pt_lim = (80,100)

    for i in range(0, len(utils.RAPIDITY_BIN), 1):
        # Retrieve the values of eta range bins (y_inf, y_sup)
        y_lim = utils.RAPIDITY_BIN[f"{i}"]

        # Cut on pt to reduce background (pt_inf, pt_sup)
        bin = 40 #binning
        mass_vs_eta(all_data,y_lim,pt_lim,bin)
        cos_vs_eta(all_data,y_lim,pt_lim,bin)

    # Compute AFB (Forward_Backward Asymmetry)
    afb(all_data, pt_lim)

    # Plot AFB for different cut in pt
    os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Z analysis')))

    # Create a txt with a list of the files from which takes values.
    pt_plot = [f"{pt_lim}"]
    for p in pt_plot:
        pathlist = Path("").glob(f"*_{p}dimuon_w.txt")
        with open(f"Afblist_{p}.txt", "w+") as outf:
            for path in pathlist:
                print(path, file=outf)

    # Dimuon Afb plots
    for p in pt_plot:
        afb_plot(f"Afblist_{p}.txt")
    os.chdir(os.path.dirname(os. getcwd()))

    # poi vedi come cambiare i programmi  per riprodurre la roba su più file

    # Elapsed time
    timer.Stop()
    print(f"Elapsed: {time.time()-start} s")
    logger.info(f"Elapsed time from the beginning is: {timer.RealTime()}(real time) seconds")
    logger.info(f"Elapsed time from the beginning is: {timer.CpuTime()}(cpu time) seconds")
