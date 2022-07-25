#include <iostream>
#include <string>
#include <ausa/json/IO.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/StringUtil.h>
#include <ausa/util/FileUtil.h>
#include <ausa/setup/DoubleSidedSiliconDetector.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Default.h>
#include <ausa/constants/Mass.h>
#include <ausa/output/OutputConvenience.h>
#include <ctime>
#include <string>
#include <TFile.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include <TFitResult.h>
#include <TF1.h>
#include <fstream>
#include <iterator>
#include <sstream>
#include <TGraphErrors.h>
#include <Math/RootFinderAlgorithms.h>
#include <Math/RootFinder.h>

using namespace::std;
using namespace ROOT::Math;

class Calibrator{
public:
    string run_number;
    TH1F *histogram;
    TH1F *sumHistogram;
    vector<double> intensity;
    vector<double> energy;
    double i_min;
    double minFittingIntensity;
    int fittingWidth;
    vector<double> peakCentroids;
    vector<double> peakCentroidSigmas;
    vector<double> guesses;
    vector<double> guessIntensities;
    vector<double> fittedPeaks;
    vector<double> energyToGuesses;
    vector<double> seedEnergies;
    vector<double> fittedIntensities;
    vector<double> tableIntensity;

    double getMinFittingIntensity() const {
        return minFittingIntensity;
    }

    void setSeedEnergies(vector<double> energies){
        seedEnergies = energies;
    }

    void setMinFittingIntensity(double minFittingIntensity) {
        Calibrator::minFittingIntensity = minFittingIntensity;
    }

    int getFittingWidth() const {
        return fittingWidth;
    }

    void setFittingWidth(int fittingWidth) {
        Calibrator::fittingWidth = fittingWidth;
    }

    void setIMin(double iMin) {
        i_min = iMin;
    }

    double getIMin() const {
        return i_min;
    }


    Calibrator(string run_number){
        this->run_number = run_number;
        i_min = 0;
        fittingWidth = 6;
        minFittingIntensity = 0;
        seedEnergies = {};
    };

    //loads histogram from the calibration file
    void loadHistogram(int clover_index){
        //open the file
        string file = "Run"+run_number+".root";
        TFile *myFile = TFile::Open(file.c_str());
        //retrieve the histogram
        string clover_toload = "h"+ to_string(clover_index)+"_Clov";
        histogram = (TH1F*)myFile->Get(clover_toload.c_str());
    }

    //find the three highest peaks in the spectrum. These can be used to create a seed calibration, which
    //is used to create the guesses for fitting the remaining peaks in the spectrum.
    vector<Int_t> findSeedPeaks(){
        vector<Int_t> peaks = {};

        //create a copy of the histogram
        TH1F *newHist = (TH1F*) histogram->Clone("newHist");

        for(int i = 0; i < seedEnergies.size(); i++){
            //find the highest peak in the spectrum
            peaks.push_back(newHist -> GetMaximumBin());
            //remove samples within +-25 samples of the highest peak.
            for(int j = peaks[i]-100; j < peaks[i]+100; j++){
                newHist ->SetBinContent(j,0);
            }
        }
        sort(peaks.begin(), peaks.end());
        return peaks;
    }

    static double linFit(double *x, double *par){
        return par[0]*x[0]+par[1];
    }

    static double poly2Fit(double *x, double *par){
        return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
    }

    static double poly3Fit(double *x, double *par){
        return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];
    }

    //find the peak above 10000 so that it can be used in seed calibration. For Zn72/Ga72 use 2201.58keV
    vector<double> findPeaksAbove10000(){
        //create a copy of the histogram
        TH1F *newHist = (TH1F*) histogram->Clone("newHist");
        //remove first 10000 bins
        for(int j = 0; j < 10000; j++){
            newHist ->SetBinContent(j,0);
        }
        double maxBin = newHist -> GetMaximumBin();
        for(int j = maxBin-100; j < maxBin+100; j++){
            newHist ->SetBinContent(j,0);
        }
        double maxBin2 = newHist -> GetMaximumBin();
        return {maxBin,maxBin2};
    }

    //create seed calibration from the seed peaks. The three highest peaks are usually 40, 121.7817, 344.27 (Eu152). Return
    //parameters of linear fit.
    vector<double> createSeedCalibration(string filename, vector<double> energyOver10000, bool useEnergyOver10000){
        auto peaks = findSeedPeaks();
        auto seedGraph  = new TGraph();
        for(int i = 0; i < seedEnergies.size(); i++){
            seedGraph ->AddPoint(peaks[i], seedEnergies[i]);
        }
        if(useEnergyOver10000){
            auto over10000 = findPeaksAbove10000();
            seedGraph -> AddPoint(over10000[0],energyOver10000[0]);
            seedGraph -> AddPoint(over10000[1],energyOver10000[1]);
        }
        auto xmax = seedGraph ->GetPointX(seedGraph->GetN()-1)+10;
        auto xmin = seedGraph ->GetPointX(0)-10;
        TF1 *func = new TF1("fit",linFit,xmin,xmax,2);
        func->SetParameters(0,0,0,0);
        TFitResultPtr fp = seedGraph->Fit("fit","s && Q && S","",xmin,xmax);
        seedGraph->Draw("AC*");
        /*auto file = new TFile(filename.c_str(), "RECREATE");
        seedGraph -> Write();
        file -> Close();*/
        return {fp ->Parameter(0), fp ->Parameter(1)};
    }

    //use this to load reference peaks.
    void loadEnergy(string filename){
        energy = {};
        intensity = {};

        std::ifstream infile (filename);

        std::string line;
        while (std::getline(infile, line))
        {
            vector<string> row_values;
            split(line, ' ', row_values);
            energy.push_back(stof(row_values[0]));
            intensity.push_back(stof(row_values[1]));
        }
    }

    // taken from http://stackoverflow.com/a/236803/248823
    void split(const std::string &s, char delim, std::vector<std::string> &elems) {
        std::stringstream ss;
        ss.str(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
    }

    //inverse lin
    static double inverseLin(double y, vector<double> params){
        return (y-params[1])/params[0];
    }

    //inverse poly
    static double inverse2Poly(double y, vector<double> params){
        auto a = params[2];
        auto b = params[1];
        auto c = params[0]-y;
        return (-b +TMath::Sqrt(b*b-4*a*c))/(2*a);
    }

    static double inversePoly(double y, vector<double> params){
        RootFinder *k = new RootFinder();
        auto *f = new TF1("fit", poly3Fit,0,16000,params.size());
        f ->SetParameter(0,params[0]-y);
        for(int i = 1; i < params.size(); i++){
            f ->SetParameter(i,params[i]);
        }
        k->Solve(*f,0,16000,1000);
        double root = k->Root();
        return root;
    }

    //Use seed calibration to create guesses for other peaks in the calibration data. Peaks with intensity lower than
    //i_min are ignored.
    void createGuesses(vector<double> seedParameters){
        guesses = {};
        for(int i = 0; i < energy.size(); i++){
            if(intensity[i] > i_min){
                guesses.push_back(inverseLin(energy[i],seedParameters));
                guessIntensities.push_back(intensity[i]);
                energyToGuesses.push_back(energy[i]);
            }
        }
    }

    static double gaussLin(double *x, double *par){
        return par[0] + par[1]*x[0] + par[2]*TMath::Gaus(x[0],par[3],par[4]);
    }

    static double doubleGaussLin(double *x, double *par){
        return par[0] + par[1]*x[0] + par[2]*TMath::Gaus(x[0],par[3],par[4], true) + par[5]*TMath::Gaus(x[0],par[6],par[7],
                                                                                                        true);
    }

    static double xGaussLin(double *x, double *par){
        int noPeaks = par[0];
        double result = par[1] + par[2]*x[0];
        for(int i = 0; i < noPeaks; i++){
            result += par[i+3]*TMath::Gaus(x[0],par[i+4],par[i+5], true);
        }
        return result;
    }

    static vector<double> findGuesses(double centroidGuess, double fittingWidth, TH1F *histogram){
        auto centroidBin = histogram -> FindBin(centroidGuess);
        auto lowerBin = histogram -> FindBin(centroidGuess - fittingWidth);
        auto upperBin = histogram -> FindBin(centroidGuess + fittingWidth);
        double ymax = 0;
        double xmax = 0;
        for(int k = lowerBin; k < upperBin; k++){
            double y = histogram->GetBinContent(k);
            if(y > ymax){
                ymax = y;
                xmax = histogram->GetBinCenter(k);
            }
        }
        return {ymax-histogram->GetBinContent(lowerBin), xmax};
    }

    //make fits to all guessed peaks. If there is another peak within the fitting width, do these simultaneously.
    void fitAllPeaks(){
        fittedPeaks = {};
        peakCentroids = {};
        peakCentroidSigmas = {};
        for(int i = 0; i < guesses.size(); i++){
            if(guessIntensities[i] < minFittingIntensity){continue;}
            //find out how many more peaks must be fitted simultaneously
            int noPeaksToFit = 0;
            bool keepRunning = true;
            for(int j = i; keepRunning; j++){
                if(!(abs(guesses[j] - guesses[j+1]) < 2*fittingWidth)){
                    keepRunning = false;
                }
                noPeaksToFit++;
            }

            //might be disturbed by Ka40 peak
            if(abs(energyToGuesses[i] - 1460) < 10){continue;}

            //create Fitting Function
            string funcname = "fit" + to_string(i);
            TF1 *func = new TF1(funcname.c_str(),xGaussLin,guesses[i]-fittingWidth,guesses[i+noPeaksToFit-1]+fittingWidth,3*noPeaksToFit+3);
            func ->FixParameter(0,noPeaksToFit);
            func ->SetParameter(1,50);
            func ->SetParameter(2,0);

            double xmax;
            double ymax;
            double sum;
            int iters;
            double xmax1;
            double xmaxlast;
            vector<double> avgs = {};

            for(int j = 0; j < noPeaksToFit; j++){
                xmax = 0;
                ymax = 0;
                sum = 0;
                iters = 0;

                for(int k = guesses[i] - fittingWidth+fittingWidth/2; k < guesses[i] + fittingWidth-fittingWidth/2; k++){
                    double y = histogram->GetBinContent(k);
                    double x = k;
                    sum += y;
                    iters++;
                    if(y > ymax){
                        ymax = y;
                        xmax = x;
                    }
                }

                bool unbroken = true;

                for(int k = xmax; unbroken; k++){
                    if(histogram->GetBinContent(k) < sum/iters){
                        avgs.push_back(k-xmax);
                        unbroken = false;
                    }
                    if(k > xmax + 800){
                        avgs.push_back(fittingWidth);
                        unbroken = false;
                    }
                }

                if(j == 0){xmax1 = xmax;}
                if(j == noPeaksToFit-1){xmaxlast = xmax;}
                func -> SetParameter(j*3+3,5*ymax);
                func -> SetParameter(j*3+4, xmax);
                func -> SetParameter(j*3+5, 5);
            }
            TFitResultPtr fp = histogram->Fit(funcname.c_str(),"+ && Q && S","",xmax1-3*avgs[0],xmaxlast+3*avgs[noPeaksToFit-1]);

            auto *fitted = histogram ->GetFunction(funcname.c_str());

            //if fit failed, continue
            if(fp-> Error(4)-guesses[i] > 5){continue;}
            if(abs(fp-> Parameter(4)-guesses[i]) > 50){continue;}

            peakCentroids.push_back(fp -> Parameter(4));
            peakCentroidSigmas.push_back(fp-> Error(4));
            fittedPeaks.push_back(energyToGuesses[i]);


            //i must be incremented to avoid fitting same peaks multiple times
            i += noPeaksToFit -1;
        }
    }

    void saveHistogram(){
        string filename = "ECalData/fittedhistograms.root";
        auto file = new TFile(filename.c_str(), "UPDATE");
        //delete old graph
        string name1 = histogram->GetName();
        string name = name1 + ";1";
        file -> Delete(name.c_str());
        histogram -> Write();
        file -> Close();
    }

    vector<double> doCalibration(){
        auto calibrationGraph = new TGraphErrors();
        calibrationGraph ->SetName(histogram->GetName());
        int lastWritten = -1;

        for(int i = 0; i < fittedPeaks.size(); i++){
            /*if(i > 0){
                if(peakCentroids[i]-calibrationGraph ->GetPointX(lastWritten) < 1000){
                    continue;
                }
            }*/
            calibrationGraph->AddPoint(peakCentroids[i],fittedPeaks[i]);
            lastWritten++;
        }

        TF1 *func = new TF1("fit",linFit,guesses[0]-10,64000,2);
        func->SetParameters(-1,0.186,3.62514e-09,0);
        TFitResultPtr fp = calibrationGraph->Fit("fit","s && Q && S","",guesses[0],guesses[guesses.size()-1]);

        string filename3 = "test.root";
        calibrationGraph->Draw("AC*");
        auto file3 = new TFile(filename3.c_str(), "UPDATE");
        calibrationGraph -> Write();
        file3 -> Close();

        //remove outliers
        for(int i = 0; i < calibrationGraph -> GetN(); i++){
            auto ch = calibrationGraph ->GetPointX(i);
            auto en = calibrationGraph ->GetPointY(i);
            auto *pars = const_cast<double*>(fp -> GetParams());
            if( abs(linFit(&ch,pars) - en) > 10){
                calibrationGraph ->RemovePoint(i);
            }
        }

        //redo calibration
        TF1 *func2 = new TF1("fit2",linFit,guesses[0]-10,64000,2);
        func2->SetParameters(-1,0.186,3.62514e-09,0);
        TFitResultPtr fp2 = calibrationGraph->Fit("fit2","s && Q && S","",guesses[0],guesses[guesses.size()-1]);

        //remove outliers
        for(int i = 0; i < calibrationGraph -> GetN(); i++){
            auto ch = calibrationGraph ->GetPointX(i);
            auto en = calibrationGraph ->GetPointY(i);
            auto *pars = const_cast<double*>(fp -> GetParams());
            if( abs(linFit(&ch,pars) - en) > 10){
                calibrationGraph ->RemovePoint(i);
            }
        }

        //redo calibration
        TF1 *func3 = new TF1("fit3",linFit,guesses[0]-10,64000,2);
        func3->SetParameters(-1,0.186,3.62514e-09,0);
        TFitResultPtr fp3 = calibrationGraph->Fit("fit3","s && Q && S","",guesses[0],guesses[guesses.size()-1]);

        auto residualGraph = new TGraphErrors();
        string histname = histogram->GetName();
        string residualname = histname + "residuals";
        residualGraph ->SetName(residualname.c_str());

        //create residualplot
        for(int i = 0; i < calibrationGraph -> GetN(); i++){
            auto ch = calibrationGraph ->GetPointX(i);
            auto en = calibrationGraph ->GetPointY(i);
            auto *pars = const_cast<double*>(fp -> GetParams());
            auto residual =linFit(&ch,pars) - en;
            residualGraph ->AddPoint(ch,residual);
        }

        string filename = "ECalData/calibrations.root";
        calibrationGraph->Draw("AC*");
        residualGraph->Draw("AC*");
        auto file = new TFile(filename.c_str(), "UPDATE");
        //delete old graph
        string name1 = calibrationGraph->GetName();
        string name = name1 + ";1";
        file -> Delete(name.c_str());
        string name2 = residualname + ";1";
        file -> Delete(name2.c_str());
        calibrationGraph -> Write();
        residualGraph -> Write();
        file -> Close();
        return {fp3 ->Parameter(0), fp3->Parameter(1)};
    }
};

class MultipleCalibrator{
public:
    vector<double> a0s;
    vector<double> a1s;
    vector<double> a2s;
    vector<double> a3s;
    vector<double> seedEnergies;
    string reffile;
    vector<double> energiesOver10000;
    bool useEnergiesOver10000;

    void setSeedEnergies(const vector<double> &seedEnergies) {
        MultipleCalibrator::seedEnergies = seedEnergies;
    }

    void setReffile(const string &reffile) {
        MultipleCalibrator::reffile = reffile;
    }

    void setEnergiesOver10000(const vector<double> &energiesOver10000) {
        MultipleCalibrator::energiesOver10000 = energiesOver10000;
    }

    void setUseEnergiesOver10000(bool useEnergiesOver10000) {
        MultipleCalibrator::useEnergiesOver10000 = useEnergiesOver10000;
    }

    Calibrator *myCalib;
    MultipleCalibrator(string run_number){
        myCalib = new Calibrator(run_number);
    }



    void calibrateAllClovers(){
        a0s = {};
        a1s = {};
        a2s = {};
        a3s = {};
        myCalib->setSeedEnergies(seedEnergies);
        for(int i = 0; i < 16; i++){
            //if(i != 10 ) continue;
            myCalib->loadHistogram(i);
            myCalib->loadEnergy(reffile);
            myCalib->setIMin(2);
            myCalib->setFittingWidth(30);
            myCalib->createGuesses(myCalib->createSeedCalibration("seed"+ to_string(i)+".root",energiesOver10000,useEnergiesOver10000));
            myCalib->fitAllPeaks();
            myCalib->saveHistogram();
            vector<double> params = myCalib->doCalibration();
            a0s.push_back(params[0]);
            a1s.push_back(params[1]);
            a2s.push_back(0);
            a3s.push_back(0);
        }
    }

    void createCalibrationFile(){
        string saveto = "ECalData/calibration.txt";
        ofstream mytxt (saveto);
        for(int i = 0; i < 16; i++){
            mytxt << "0\t" << i << "\t" << a0s[i] << "\t" << a1s[i] << "\t" << a2s[i] << "\t" << a3s[i] << "\t" << 0 << "\n";
        }
        //also Si calibration?
        mytxt << "4\t" << 0 << "\t" << 233.687 << "\t" << 2.20288 << "\t" << 7.21896E-05 << "\t" << 0 << "\t" << 0 << "\n";
        mytxt.close();
    }

};

int main() {
    /*auto myMulti = new MultipleCalibrator("088Old");
    myMulti ->setSeedEnergies({144.6,191.5,629.96,834.1});
    myMulti ->setUseEnergiesOver10000(true);
    myMulti ->setEnergiesOver10000({2201.58,2507.75});
    myMulti ->setReffile("refsGa.txt");
    myMulti -> calibrateAllClovers();
    myMulti -> createCalibrationFile();*/
    auto myMulti = new MultipleCalibrator("091");
    myMulti ->setSeedEnergies({39.9047,121.7817,344.2789});
    myMulti ->setUseEnergiesOver10000(false);
    myMulti ->setEnergiesOver10000({2201.58,2507.75});
    myMulti ->setReffile("refsEu.txt");
    myMulti -> calibrateAllClovers();
    myMulti -> createCalibrationFile();
}
