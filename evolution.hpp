//
// Created by polya on 4/2/24.
//

#ifndef OS_DCE_CPP_EVOLUTION_HPP
#define OS_DCE_CPP_EVOLUTION_HPP
#include <string>
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <regex>
#include <array>
#include <boost/filesystem.hpp>
#include <armadillo>


namespace fs = boost::filesystem;
using namespace std::complex_literals;
const auto PI=std::numbers::pi;

//This subroutine computes evolution using operator splitting and particle number
class os_DCE_Evolution {


public:
    /// This constructor initializes all parameters
    /// @param group group number of parameters
    /// @param row row number of parameters
    os_DCE_Evolution(const int &group, const int &row) {
        //
        this->groupNum = group;
        this->rowNum = row;

        this->parseCSV(group, row);
        for (int n1 =0;n1<N1;n1++){
            this->x1ValsAll.push_back(-L1+dx1*n1);
        }
        for (int n2=0;n2<N2;n2++){
            this->x2ValsAll.push_back(-L2+dx2*n2);
        }

        for(const auto& val: x1ValsAll){
            x1ValsAllSquared.push_back(std::pow(val,2));
        }
        for(const auto &val:x2ValsAll){
            x2ValsAllSquared.push_back(std::pow(val,2));
        }

        for(int n1=0;n1<static_cast<int>(N1/2);n1++){
            k1ValsAll.push_back(2*PI*static_cast<double >(n1)/(2.0*L1));
        }
        for(int n1=static_cast<int>(N1/2);n1<N1;n1++){
            k1ValsAll.push_back(2*PI*static_cast<double >(n1-N1)/(2.0*L1));
        }

        for(const auto&val: k1ValsAll){
            k1ValsAllSquared.push_back(std::pow(val,2));
        }

        for(int n2=0;n2<static_cast<int>(N2/2);n2++){
            k2ValsAll.push_back(2*PI*static_cast<double >(n2)/(2.0*L2));
        }
        for(int n2=static_cast<int >(N2/2);n2<N2;n2++){
            k2ValsAll.push_back(2*PI*static_cast<double >(n2-N2)/(2.0*L2));
        }
        for(const auto &val:k2ValsAll){
            k2ValsAllSquared.push_back(std::pow(val,2));
        }


        for(int fls=0;fls<flushNum;fls++){
            int startingInd=fls*stepsPerFlush;
            for(int j=0;j<stepsPerFlush;j++){
                double indTmp=static_cast<double >(startingInd+j);
                this->timeValsAll.push_back(indTmp*dt);
            }
        }
//        printVec(timeValsAll);
//        std::cout<<"dt="<<dt<<std::endl;
//        std::cout<<"dx1="<<dx1<<std::endl;
//        std::cout<<"dx2="<<dx2<<std::endl;
//        std::cout<<"M="<<M<<std::endl;
//        std::cout<<"x1vec has size "<<x1ValsAll.size()<<std::endl;
//        std::cout<<"x2vec has size "<<x2ValsAll.size()<<std::endl;

    }//end of constructor





public:
    int jH1 = -1;
    int jH2 = -1;
    double g0 = 0;
    double omegam = 0;
    double omegap=0;
    double omegac = 0;
    double er = 0;
    double thetaCoef = 0;
    int groupNum = -1;
    int rowNum = -1;
    double theta=0;
    double lmd=0;
    double Deltam=0;

    int N1=10;//000;
    int N2=12;//000;

    double L1=0.5;
    double L2=0.8;
    double dx1=0;
    double dx2=0;

    double dtEst=0.0001;
    double tFlushStart=0;
    double tFlushStop=0.001;
    double tTotPerFlush=tFlushStop-tFlushStart;
    int flushNum=10;
    int stepsPerFlush=static_cast<int>(std::ceil(tTotPerFlush/dtEst));
    double dt=tTotPerFlush/static_cast<double >(stepsPerFlush);
    std::vector<double> timeValsAll;


    std::vector<double> x1ValsAll;
    std::vector<double> x2ValsAll;
    std::vector<double> k1ValsAll;
    std::vector<double> k2ValsAll;
    std::vector<double> x1ValsAllSquared;
    std::vector<double> x2ValsAllSquared;
    std::vector<double> k1ValsAllSquared;
    std::vector<double> k2ValsAllSquared;

//    arma::cx_dvec psi0;
    arma::cx_dmat psi0;



public:

    /// @param group group number
    /// @param row row number
    ///parse csv with group number group
    void parseCSV(const int &group, const int &row);

    ///
    /// @param cmd python execution string
    /// @return signal from the python
    static std::string execPython(const char *cmd);

    template<class T>
    static void printVec(const std::vector <T> &vec) {
        for (int i = 0; i < vec.size() - 1; i++) {
            std::cout << vec[i] << ",";
        }
        std::cout << vec[vec.size() - 1] << std::endl;
    }

    ///initialize wavefunction serially
    void initPsiSerial();



    ///
    /// @param n1 index of x1
    /// @return wavefunction of photon at n1
    double f1(int n1);

    ///
    /// @param n2 index of x2
    /// @return wavefunction of phonon at n2
    double f2(int n2);

    ///
    /// @param n1 index for x1n1
    /// @param t time
    /// @return coefficient for evolution using H3
    double f(int n1, double t);

    ///
    /// @param j time step
    /// @param psi wavefunction at the beginning of the time step j
    /// @return
    arma::cx_dmat evolution1Step(const int&j, const arma::cx_dmat& psi);

};


#endif //OS_DCE_CPP_EVOLUTION_HPP
