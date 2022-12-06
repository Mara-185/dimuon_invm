# TEST
import ROOT
import unittest
import math
import numpy as np
import utils
import dimuon_ana
#from utils import *
#from dimuon_ana import *

"""This script generate some events of mumu in Z in order to test the behaviour
    of the main functions defined in \"dimuon_ana.py\".
    The Z boson is chosed to test also the functions used to study the weak
    mixing angle from its angular properties. """

# Generate a root file for the testing

#############################
# TEST CONSIDERANDO UNA CLASSE PARTICLE DA CUI EREDITANO ALTRI ESEMPI
# POI RDATAFRAME ...COME?
# POI LEPTONS ANALISYS E BLA BLA????

def create_dataset(N):

    # Enable parallel analysis
    nthreads = 4
    ROOT.ROOT.EnableImplicitMT(nthreads)

    # Create the RDataFrame of lenght N
    rdf = ROOT.RDataFrame(N)
    ROOT.gRandom.SetSeed(1)
    nbin = math.floor(math.sqrt(N))

    # By JITting define two functions:
    # - "pseudorapidity": which generate a sample from two different gaussian,
    # centered in two different values, obtained by fitting the Pseudorapidity
    # distribution from real data;
    # - "mass": which generate the mass distribution from a gaussian peaked in
    # 90.8 and from an exponential distribution to simulate the background.
    # The initial parameters are obtained from the fit of real values, too.

    ROOT.gInterpreter.Declare("""
    double pseudorapidity(){

        double a;
        a = gRandom->Uniform(0,1);
        if(a<0.5){
            return gRandom->Gaus(-2.47, 1.41);}
        else{
            return gRandom->Gaus(2.47, 1.41);}
    }

    double mass(){
        double b;
        b = gRandom->Uniform(0,1);
        if(b<0.57){
            return gRandom->Gaus(90.8, 2.7);}
        else{
            return (gRandom->Exp(33.3))+60;}
    }

    """)

    # Define the column of the dataset
    rdf_m = rdf.Define("Dimuon_mass", "mass()").\
        Define("Dimuon_phi", "gRandom->Uniform(-3.14, 3.14)").\
        Define("Dimuon_pt", "gRandom->Landau(7.14, 3.88)").\
        Define("Dimuon_eta", "pseudorapidity()")

    del rdf

    #Verify the distributions
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("mass", "mass",
        nbin, 0, 110), "Dimuon_mass")
    h1 = rdf_m.Histo1D(ROOT.RDF.TH1DModel("phi", "phi",
        nbin, -3.5, 3.5), "Dimuon_phi")
    h2 = rdf_m.Histo1D(ROOT.RDF.TH1DModel("pt", "pt",
        nbin, 0, 120), "Dimuon_pt")
    h3 = rdf_m.Histo1D(ROOT.RDF.TH1DModel("eta", "eta",
        nbin, -10, 10), "Dimuon_eta")

    branchlist = ["Dimuon_mass", "Dimuon_phi", "Dimuon_pt", "Dimuon_eta"]
    data_cached = rdf_m.Snapshot("data", "data.root", branchlist).Cache()

    c = ROOT.TCanvas()
    h.Draw()
    c.SaveAs("test_mass.png")
    h1.Draw()
    c.SaveAs("test_phi.png")
    h2.Draw()
    c.SaveAs("test_pt.png")
    h3.Draw()
    c.SaveAs("test_eta.png")

    return data_cached


def create_dataset2(N):

    # Enable parallel analysis
    nthreads = 4
    ROOT.ROOT.EnableImplicitMT(nthreads)

    # Create the RDataFrame of lenght N
    rdf = ROOT.RDataFrame(N)
    ROOT.gRandom.SetSeed(1)
    nbin = math.floor(math.sqrt(N))


    ROOT.gInterpreter.Declare("""

    ROOT::VecOps::RVec<double> phi2(){
        double a,phi0, phi1;
        phi0 = gRandom->Uniform(-3.14,3.14);
        phi1 = gRandom->Uniform(-3.14, 3.14);

        ROOT::VecOps::RVec<double> Phi{phi0,phi1};
        return Phi;
    }

    ROOT::VecOps::RVec<double> mass2(){
        double a,mass0, mass1;
        mass0 = gRandom->Uniform(0.1056583, 0.105659);
        mass1 = gRandom->Uniform(0.1056583, 0.105659);

        ROOT::VecOps::RVec<double> Mass{mass0,mass1};
        return Mass;
    }


    ROOT::VecOps::RVec<double> pt2(){
        double a,b, pt0, pt1;
        a = gRandom->Uniform(0,1);
        if((a>=0.7)&&(a<0.85)){
            pt0= gRandom->Landau(0, .8);}
        else {
            if(a<0.4){
                pt0= gRandom->Landau(10, 1);}
            else {
                if((a>=0.4)&&(a<0.7)){
                    pt0= gRandom->Landau(15, 1);}
                else {
                    if((a>=0.85)&&(a<1)){
                        pt0= gRandom->Landau(43, 4);}}}}

        b = gRandom->Uniform(0,1);
        if((b>=0.7)&&(b<0.85)){
            pt1= gRandom->Landau(0, .8);}
        else {
            if(b<0.4){
                pt1= gRandom->Landau(10, 1);}
            else {
                if((b>=0.4)&&(b<0.7)){
                    pt1= gRandom->Landau(15, 1);}
                else {
                    if((b>=0.85)&&(b<1)){
                        pt1= gRandom->Landau(43, 4);}}}}
        ROOT::VecOps::RVec<double> Pt{pt0,pt1};
        return Pt;
    }

    ROOT::VecOps::RVec<double> eta2(){

        double a, b, eta0, eta1;
        a = gRandom->Uniform(0,1);
        if((a>=0.03)&&(a<.35)){
            eta0= gRandom->Uniform(-2.5, 2.5);}
        else {
            if(a<0.03){
                eta0= gRandom->Uniform(-10,10);}
            else {
                if((a>=0.35)&&(a<0.65)){
                    eta0= gRandom->Uniform(-1.8,1.8);}
                else {
                    if((a>=.65)&&(a<.92)){
                        eta0= gRandom->Uniform(-1.18, 1.18);}
                    else {
                        if((a>=.92)&&(a<0.95)){
                            eta0= gRandom->Uniform(0.3,0.8);}
                        else {
                            if((a>=.95)&&(a<.98)){
                                eta0= gRandom->Uniform(-0.3,-0.8);}
                            else {
                                if((a>=.98)&&(a<.99)){
                                    eta0= gRandom->Uniform(-0.2,-0.01);}
                                else {
                                    if((a>=.99)&&(a<1)){
                                        eta0= gRandom->Uniform(0.01,0.2);}}}}}}}}

        b = gRandom->Uniform(0,1);
        if((b>=0.03)&&(b<.35)){
            eta1= gRandom->Uniform(-2.5, 2.5);}
        else {
            if(b<0.03){
                eta1= gRandom->Uniform(-10,10);}
            else {
                if((b>=0.35)&&(b<0.65)){
                    eta1= gRandom->Uniform(-1.8,1.8);}
                else {
                    if((b>=.65)&&(b<.92)){
                        eta1= gRandom->Uniform(-1.18, 1.18);}
                    else {
                        if((b>=.92)&&(b<0.95)){
                            eta1= gRandom->Uniform(0.3,0.8);}
                        else {
                            if((b>=.95)&&(b<.98)){
                                eta1= gRandom->Uniform(-0.3,-0.8);}
                            else {
                                if((b>=.98)&&(b<.99)){
                                    eta1= gRandom->Uniform(-0.2,-0.01);}
                                else {
                                    if((b>=.99)&&(b<1)){
                                        eta1= gRandom->Uniform(0.01,0.2);}}}}}}}}
        ROOT::VecOps::RVec<double> Eta{eta0,eta1};
        return Eta;
    }

    """)

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


    # Define the column of the dataset
    rdf_m = rdf.Define("Muon_mass", "mass2()").\
        Define("Muon_phi", "phi2()").\
        Define("Muon_pt", "pt2()").\
        Define("Muon_eta", "eta2()")

    rdf_dimu = rdf_m.\
        Define("Dimuon_mass", \
            "dilepton_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
            Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[3]").\
        Define("Dimuon_pt",\
            "dilepton_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[0]").\
        Define("Dimuon_eta", \
            "dilepton_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[1]").\
        Define("Dimuon_phi",\
            "dilepton_vec(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[2]").\
        Define("Dimuon_cos",\
            "cos_rapidity(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[0]").\
        Define("Dimuon_y",\
            "cos_rapidity(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], \
             Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])[1]")


    #Verify the distributions
    h = rdf_m.Histo1D(ROOT.RDF.TH1DModel("mass", "mass",
        nbin, 0.1, 0.2), "Muon_mass")
    h1 = rdf_m.Histo1D(ROOT.RDF.TH1DModel("phi", "phi",
        nbin, -3.5, 3.5), "Muon_phi")
    h2 = rdf_m.Histo1D(ROOT.RDF.TH1DModel("pt", "pt",
        nbin, 0, 120), "Muon_pt")
    h3 = rdf_m.Histo1D(ROOT.RDF.TH1DModel("eta", "eta",
        nbin, -10, 10), "Muon_eta")

    branchlist = ["Muon_mass", "Muon_phi", "Muon_pt", "Muon_eta"]#, "Muon_cos", "Muon_y"]
    data_cached = rdf_m.Snapshot("Events", "events.root", branchlist)#.Cache()
    branchlist_mu = ["Dimuon_mass", "Dimuon_pt", "Dimuon_eta", "Dimuon_phi", \
        "Dimuon_cos", "Dimuon_y"]
    rdf_dimu.Snapshot("dimuon", "dimuon.root", branchlist_mu)

    c = ROOT.TCanvas()
    h.Draw()
    c.SaveAs("test_mass2.png")
    h1.Draw()
    c.SaveAs("test_phi2.png")
    h2.Draw()
    c.SaveAs("test_pt2.png")
    h3.Draw()
    c.SaveAs("test_eta2.png")


# Actual test

# 1. dilepton_ana
# 2. fit
# 3. PROPERTIES
# 4. Z analysis

# class TestDimuon(unittest.TestCase):
#
#     def test_fit(data):
#         resonance_fit(data, "Z")
#
#     def test_prop(data, cached):
#         resonance_prop(data, cached , "Z")
#
#     def test_leptons(data):
#         leptons_analysis(data)


def test_fit(data):
    resonance_fit(data, "Z")

def test_prop(data, cached):
    resonance_prop(data, cached , "Z")

def test_leptons(data):
    leptons_analysis(data)


def test_weight(data):
    weight(data)

def test_afb(infile):
    afb(infile)


if __name__=="__main__":

    # Start the the timer
    timer = ROOT.TStopwatch()
    timer.Start()

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

    # Creating the logger and setting its level
    utils.set_logger("Dimuon_invm_tests")

    # Enable parallel analysis
    nthreads = 4
    ROOT.ROOT.EnableImplicitMT(nthreads)

    logger.info("Starting the tests...")

    # Define number of events to generate
    N = 100000
    outfile = "data.root"
    #data = create_dataset(N)
    logger.info("The dataset is created.")

    # Begin the tests
    logger.info("Start test on fit function.")

    # Fit
    os.makedirs("Fit", exist_ok=True)
    logger.debug("The new directory \"Fit\" is created")
    #test_fit(outfile)

    # Properties
    os.makedirs("Properties", exist_ok=True)
    logger.debug("The new directory \"Properties\" is created")
    #test_prop(outfile, data)

##############################################
    #unittest.main()


    # Z ANALYSIS
    os.makedirs("Z analysis", exist_ok=True)
    logger.debug("The new directory \"Z analysis\" is created")
    # N2 = 100000
    # outfile = "dimuon.root"
    # outfile_w = "dimuon_w.root"
    # create_dataset2(N2)
    #
    # test_weight(outfile)
    # test_afb(outfile_w)
    #
    # # Plot AFB for different cut in pt
    #
    # os.chdir(os.path.abspath(os.path.join(os.sep,f'{os.getcwd()}', 'Z analysis')))
    # pt_plot = ["(10,120)"]
    # for p in pt_plot:
    #     pathlist = Path("").glob(f"*_{p}dimuon_w.txt")
    #     with open(f"Afblist_{p}.txt", "w+") as outf:
    #         for path in pathlist:
    #             print(path, file=outf)
    #
    # ##Dimuon Afb plots
    # for p in pt_plot:
    #     afb_plot(f"Afblist_{p}.txt")
    # os.chdir(os.path.dirname(os. getcwd()))




    timer.Stop()
    logger.info(f"Elapsed time from the beginning is: {N} -> {timer.RealTime()}(real time) seconds")
