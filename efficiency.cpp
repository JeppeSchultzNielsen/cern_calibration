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

using namespace::std;

class EfficiencyCalibrator{
public:
    string run_number;
    TH1F *histogram;
    TH1F *sumHistogram;
    vector<double> tableIntensity;

    void setIMin(double iMin) {
        i_min = iMin;
    }

    vector<double> tableEnergy;
    double i_min;
    double minFittingIntensity;
    int fittingWidth;

    void setFittingWidth(int fittingWidth) {
        EfficiencyCalibrator::fittingWidth = fittingWidth;
    }

    vector<double> tableIntensityAboveImin;
    vector<double> tableEnergyAboveImin;
    vector<double> fittedIntensity;
    vector<double> fittedTableIntensity;
    vector<double> fittedTableEnergy;

    vector<vector<double>> individualFittedIntensity;
    vector<vector<double>> individualFittedTableIntensity;
    vector<vector<double>> individualFittedTableEnergy;


    EfficiencyCalibrator(string run_number){
        this->run_number = run_number;
        i_min = 0;
        fittingWidth = 6;
        minFittingIntensity = 0;
    };

    //assume calibrated clovers. Combine all histograms
    void combineAllHistograms(){
        //open the file
        string file = "Run"+run_number+".root";
        unique_ptr<TFile> myFile(TFile::Open(file.c_str()));
        //retrieve the first histogram
        string clover_toload = "h"+ to_string(0)+"_Clov";
        sumHistogram = (TH1F*)myFile->Get(clover_toload.c_str());
        //add all other histograms into sumHistogram
        for(int i = 0; i < 16; i++){
            string clover_toload = "h"+ to_string(i)+"_Clov";
            sumHistogram -> Add((TH1F*)myFile->Get(clover_toload.c_str()));
        }
        string filename = "combinedHistogram.root";
        auto file2 = new TFile(filename.c_str(), "RECREATE");
        file2 -> cd();
        sumHistogram -> Write();
        file2 -> Close();
        myFile -> Close();
    }

    void createHistogram(double lower, double upper, int bins){
        string file = "Run"+run_number+".root";
        unique_ptr<TFile> myFile(TFile::Open(file.c_str()));
        TTree *t = (TTree *) myFile->Get("ids");
        auto entries = t -> GetEntries();
        double_t Clov_En[16];
        t->SetBranchAddress("Energy_Clov", &Clov_En);

        double currentEn;
        sumHistogram = new TH1F("histogram","histogram",bins,lower,upper);
        for(int i = 0; i < entries; i++){
            t ->GetEntry(i);
            for(int j = 0; j < 16; j++){
                currentEn = Clov_En[j];
                if(currentEn > 0.1){
                    sumHistogram -> Fill(currentEn);
                }
            }
        }
        string filename = "EfficiencyCalData/combinedHistogram.root";
        auto file2 = new TFile(filename.c_str(), "RECREATE");
        file2 -> cd();
        sumHistogram -> Write();
        file2 -> Close();
        myFile -> Close();
    }

    void createHistogramsIndiviually(double lower, double upper, int bins){
        string file = "Run"+run_number+".root";
        unique_ptr<TFile> myFile(TFile::Open(file.c_str()));
        TTree *t = (TTree *) myFile->Get("ids");
        auto entries = t -> GetEntries();
        double_t Clov_En[16];
        t->SetBranchAddress("Energy_Clov", &Clov_En);

        double currentEn;
        TH1F *histograms[16];
        for(int j = 0; j < 16; j++){
            string histname = "h"+to_string(j);
            histograms[j] = new TH1F(histname.c_str(),histname.c_str(),bins,lower,upper);
        }

        for(int i = 0; i < entries; i++){
            t ->GetEntry(i);
            for(int j = 0; j < 16; j++){
                currentEn = Clov_En[j];
                if(currentEn > 0.1){
                    histograms[j] -> Fill(currentEn);
                }
            }
        }
        string filename = "EfficiencyCalData/individualHistograms"+run_number+".root";
        auto file2 = new TFile(filename.c_str(), "RECREATE");
        file2 -> cd();
        for(int j = 0; j < 16; j++){
            histograms[j] -> Write();
        }
        file2 -> Close();
        myFile -> Close();
    }

    void loadEnergy(string filename){
        tableEnergy = {};
        tableIntensity = {};

        std::ifstream infile (filename);

        std::string line;
        while (std::getline(infile, line))
        {
            vector<string> row_values;
            split(line, ' ', row_values);
            tableEnergy.push_back(stof(row_values[0]));
            tableIntensity.push_back(stof(row_values[1]));
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

    void findPeaksOverImin(){
        for(int i = 0; i < tableEnergy.size(); i++) {
            if (tableIntensity[i] > i_min) {
                tableEnergyAboveImin.push_back(tableEnergy[i]);
                auto intensity = tableIntensity[i];
                for (int j = 0; j < tableEnergy.size(); j++) {
                    if (j == i) { continue; }
                    auto difference = abs(tableEnergy[i] - tableEnergy[j]);
                    if (difference < 1) {
                        intensity += tableIntensity[j];
                    }
                    if (tableEnergy[j] > tableEnergy[i] + 10) {
                        j = tableEnergy.size();
                    }
                }
                tableIntensityAboveImin.push_back(intensity);
            }
        }
    }

    static double xGaussLin(double *x, double *par){
        int noPeaks = par[0];
        double result = par[1] + par[2]*x[0];
        for(int i = 0; i < noPeaks; i++){
            result += par[i+3]*TMath::Gaus(x[0],par[i+4],par[i+5], true);
        }
        return result;
    }

    static double gaussLin(double *x, double *par){
        return par[0] + par[1]*x[0] + par[2]*TMath::Gaus(x[0],par[3],par[4], true);
    }

    void fitAllPeaksEfficiency(){
        string file = "EfficiencyCalData/combinedHistogram.root";
        unique_ptr<TFile> myFile(TFile::Open(file.c_str()));
        string clover_toload = "histogram";
        sumHistogram = (TH1F*)myFile->Get(clover_toload.c_str());
        fittedIntensity = {};
        fittedTableEnergy = {};
        fittedTableIntensity = {};

        for(int i = 0; i < tableEnergyAboveImin.size(); i++){
            //find out how many more peaks must be fitted simultaneously
            int noPeaksToFit = 0;
            bool keepRunning = true;
            for(int j = i; keepRunning; j++){
                if(!(abs(tableEnergyAboveImin[j] - tableEnergyAboveImin[j+1]) < 2*fittingWidth)){
                    keepRunning = false;
                }
                noPeaksToFit++;
            }
            if(noPeaksToFit > 1){i += noPeaksToFit -1; continue;}
            TF1 *func = new TF1("fit",gaussLin,tableEnergyAboveImin[i]-fittingWidth,tableEnergyAboveImin[i+noPeaksToFit-1]+fittingWidth,5);
            auto params = findGuesses(tableEnergyAboveImin[i], fittingWidth, sumHistogram);
            func ->SetParameters(0,0,params[0],params[1],fittingWidth/5);
            TFitResultPtr fp = sumHistogram->Fit("fit","+ && Q && S","",tableEnergyAboveImin[i]-fittingWidth,tableEnergyAboveImin[i+noPeaksToFit-1]+fittingWidth);

            fittedIntensity.push_back(fp -> Parameter(2));
            fittedTableIntensity.push_back(tableIntensityAboveImin[i]);
            fittedTableEnergy.push_back(tableEnergyAboveImin[i]);

            //i must be incremented to avoid fitting same peaks multiple times
            i += noPeaksToFit -1;
        }
        string file2 = "EfficiencyCalData/combinedHistogramFit.root";
        auto myFile2 = new TFile(file2.c_str(), "RECREATE");
        myFile2 -> cd();
        sumHistogram -> Write();
        myFile2->Close();
        myFile -> Close();
    }

    void saveData(string filename){
        string saveto = filename;
        ofstream mytxt (saveto);
        for(int i = 0; i < fittedIntensity.size(); i++){
            mytxt << fittedTableEnergy[i] << "\t" << fittedTableIntensity[i] << "\t" << fittedIntensity[i] << "\n";
        }
        mytxt.close();
    }

    void fitAllPeaksIndividuallyEfficiency(){
        individualFittedIntensity = {};
        individualFittedTableEnergy = {};
        individualFittedTableIntensity = {};

        string file = "EfficiencyCalData/individualHistograms"+run_number+".root";
        unique_ptr<TFile> myFile(TFile::Open(file.c_str()));

        for(int k = 0; k < 16; k++){
            string clover_toload = "h"+to_string(k);
            histogram = (TH1F *) myFile->Get(clover_toload.c_str());
            vector<double> fittedIntensity1 = {};
            vector<double> fittedTableEnergy1 = {};
            vector<double> fittedTableIntensity1 = {};

            for (int i = 0; i < tableEnergyAboveImin.size(); i++) {
                //find out how many more peaks must be fitted simultaneously
                int noPeaksToFit = 0;
                bool keepRunning = true;
                for (int j = i; keepRunning; j++) {
                    if (!(abs(tableEnergyAboveImin[j] - tableEnergyAboveImin[j + 1]) < 2 * fittingWidth)) {
                        keepRunning = false;
                    }
                    noPeaksToFit++;
                }
                if (noPeaksToFit > 1) {
                    i += noPeaksToFit - 1;
                    continue;
                }
                TF1 *func = new TF1("fit", gaussLin, tableEnergyAboveImin[i] - fittingWidth,
                                    tableEnergyAboveImin[i + noPeaksToFit - 1] + fittingWidth, 5);
                auto params = findGuesses(tableEnergyAboveImin[i], fittingWidth, histogram);
                func->SetParameters(0, 0, 2*params[0], params[1], fittingWidth / 10);

                TFitResultPtr fp = histogram->Fit("fit", "+ && Q && S", "", params[1] - fittingWidth,
                                                  params[1] + fittingWidth);

                fittedIntensity1.push_back(fp->Parameter(2));
                fittedTableIntensity1.push_back(tableIntensityAboveImin[i]);
                fittedTableEnergy1.push_back(tableEnergyAboveImin[i]);

                //i must be incremented to avoid fitting same peaks multiple times
                i += noPeaksToFit - 1;
            }

            individualFittedIntensity.push_back(fittedIntensity1);
            individualFittedTableIntensity.push_back(fittedTableIntensity1);
            individualFittedTableEnergy.push_back(fittedTableEnergy1);

            string file2 = "EfficiencyCalData/individualFits"+run_number+".root";
            auto myFile2 = new TFile(file2.c_str(), "UPDATE");
            //delete old graph
            string name1 = histogram->GetName();
            string name = name1 + ";1";
            myFile2 -> Delete(name.c_str());
            myFile2->cd();
            histogram->Write();
            myFile2->Close();
        }
        myFile->Close();
    }

    void saveIndividualEfficiency(string filename){
        string saveto = filename;
        ofstream mytxt (saveto);
        for(int i = 0; i < individualFittedIntensity[0].size(); i++){
            double peakSum = 0;
            for(int j = 0; j < 16; j++){
                peakSum += individualFittedIntensity[j][i];
            }
            mytxt << individualFittedTableEnergy[0][i] << "\t" << individualFittedTableIntensity[0][i] << "\t" << peakSum << "\n";
        }
        mytxt.close();
    }

    void printTimeDiff(){
        string file = "Run"+run_number+".root";
        unique_ptr<TFile> myFile(TFile::Open(file.c_str()));
        TTree *t = (TTree *) myFile->Get("ids");
        auto entries = t -> GetEntries();
        ULong64_t timestamp;
        t->SetBranchAddress("Timestamp", &timestamp);
        t -> GetEntry(1);
        auto t_1 = timestamp;
        t -> GetEntry(entries - 1);
        auto t_2 = timestamp;
        cout << t_2 - t_1 << endl;
    }
};

int main() {
    auto myCalib = new EfficiencyCalibrator("091");
    //myCalib ->createHistogramsIndiviually(0,3000,12000);
    myCalib ->setIMin(1);
    myCalib ->setFittingWidth(10);
    myCalib ->loadEnergy("refsEu.txt");
    myCalib -> findPeaksOverImin();
    myCalib -> fitAllPeaksIndividuallyEfficiency();
    myCalib ->saveIndividualEfficiency("EfficiencyCalData/individualEu.txt");

    /*myCalib ->setIMin(1);
    myCalib ->setFittingWidth(10);
    myCalib ->loadEnergy("refs.txt");
    myCalib -> findPeaksOverImin();
    myCalib ->createHistogram(0,3000,12000);
    myCalib -> fitAllPeaksEfficiency();
    myCalib -> saveData("efficiencyCalibration.txt");
    myCalib -> printTimeDiff();*/
}