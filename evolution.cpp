//
// Created by polya on 4/2/24.
//

#include "evolution.hpp"

/// @param group group number
/// @param row row number
///parse csv with group number group
void os_DCE_Evolution::parseCSV(const int &group, const int &row){
    std::string commandToReadCSV="python3 readCSV.py "+std::to_string(group)+" "+std::to_string(row);

    std::string result=this->execPython(commandToReadCSV.c_str());
//    std::cout<<result<<std::endl;

    std::regex pattern_j1H("j1H(\\d+)j2H");
    std::smatch  match_j1H;
    if (std::regex_search(result,match_j1H,pattern_j1H)){
        this->jH1=std::stoi(match_j1H[1].str());
    }

    std::regex pattern_j2H("j2H(\\d+)g0");
    std::smatch match_j2H;
    if (std::regex_search(result,match_j2H,pattern_j2H)){
        this->jH2=std::stoi(match_j2H[1].str());
    }

    std::regex pattern_g0("g0([+-]?\\d+(\\.\\d+)?)omegam");
    std::smatch match_g0;
    if (std::regex_search(result,match_g0,pattern_g0)){
        this->g0=std::stod(match_g0[1].str());
    }


    std::regex pattern_omegam("omegam([+-]?\\d+(\\.\\d+)?)omegap");
    std::smatch match_omegam;
    if (std::regex_search(result,match_omegam,pattern_omegam)){
        this->omegam=std::stod(match_omegam[1].str());
    }


    std::regex pattern_omegap("omegap([+-]?\\d+(\\.\\d+)?)omegac");
    std::smatch match_omegap;
    if (std::regex_search(result,match_omegap,pattern_omegap)){
        this->omegap=std::stod(match_omegap[1].str());
    }

    std::regex pattern_omegac("omegac([+-]?\\d+(\\.\\d+)?)er");
    std::smatch match_omegac;
    if (std::regex_search(result,match_omegac,pattern_omegac)){
        this->omegac=std::stod(match_omegac[1].str());
    }

    std::regex pattern_er("er([+-]?\\d+(\\.\\d+)?)thetaCoef");
    std::smatch match_er;
    if(std::regex_search(result,match_er,pattern_er)){
        this->er=std::stod(match_er[1].str());
    }

    std::regex pattern_thetaCoef("thetaCoef([+-]?\\d+(\\.\\d+)?)");
    std::smatch  match_thetaCoef;
    if (std::regex_search(result,match_thetaCoef,pattern_thetaCoef)){
        this->thetaCoef=std::stod(match_thetaCoef[1].str());
    }
//    std::cout<<"jH1="<<jH1<<std::endl;
//
//    std::cout<<"jH2="<<jH2<<std::endl;
//
//    std::cout<<"g0="<<g0<<std::endl;
//
//    std::cout<<"omegam="<<omegam<<std::endl;
//
//    std::cout<<"omegac="<<omegac<<std::endl;
//
//    std::cout<<"omegap="<<omegap<<std::endl;
//
//    std::cout<<"er="<<er<<std::endl;
//
//    std::cout<<"thetaCoef="<<thetaCoef<<std::endl;
    double e2r=std::pow(er,2);
    double eM2r=1/e2r;
    this->Deltam=this->omegam-this->omegap;
    this->lmd=(e2r-eM2r)/(e2r+eM2r)*Deltam;
    this->theta=thetaCoef*PI;
//      std::cout<<"lambda="<<lmd<<std::endl;
//      std::cout<<"theta="<<theta<<std::endl;
//      std::cout<<"Deltam="<<Deltam<<std::endl;
//    double height1=0.5;
//
//    double width1=std::sqrt(-2.0*std::log(height1)/omegac);
//    double minGrid1=width1/20.0;
//    this->N1=static_cast<int>(std::ceil(L1*2/minGrid1));
//    if(N1%2==1){
//        N1+=1;//make sure N1 is even
//    }

    std::cout<<"N1="<<N1<<std::endl;
    dx1=2*L1/(static_cast<double>(N1));
    dx2=2*L2/(static_cast<double >(N2));
//    std::cout<<"dt="<<dt<<std::endl;
}

///
/// @param cmd python execution string
/// @return signal from the python
 std::string os_DCE_Evolution::execPython(const char *cmd){
    std::array<char, 4096> buffer; // Buffer to store command output
    std::string result; // String to accumulate output

    // Open a pipe to read the output of the executed command
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }

    // Read the output a chunk at a time and append it to the result string
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }

    return result; // Return the accumulated output


}




///
/// @param n1 index of x1
/// @return wavefunction of photon at n1
double os_DCE_Evolution::f1(int n1){
    double x1TmpSquared=x1ValsAllSquared[n1];
    double x1Tmp=x1ValsAll[n1];

    double valTmp = std::exp(-0.5 * omegac * x1TmpSquared)
                    * std::hermite(this->jH1, std::sqrt(omegac) * x1Tmp);


    return valTmp;


}


///
/// @param n2 index of x2
/// @return wavefunction of phonon at n2
double os_DCE_Evolution::f2(int n2){
    double x2TmpSquared=x2ValsAllSquared[n2];
    double x2Tmp=x2ValsAll[n2];

    double valTmp=std::exp(-0.5 * omegam * x2TmpSquared)
                  *std::hermite(this->jH2,std::sqrt(omegam)*x2Tmp);

    return valTmp;

}


///initialize wavefunction serially
void os_DCE_Evolution::initPsiSerial(){

    arma::cx_dcolvec vec1(N1);
    arma::cx_drowvec vec2(N2);
    for(int n1=0;n1<N1;n1++){
        vec1(n1)= f1(n1);
    }
    for(int n2=0;n2<N2;n2++){
        vec2(n2)= f2(n2);
    }
    this->psi0=arma::kron(vec1,vec2);
    this->psi0/=arma::norm(psi0,2);

//    std::cout<<psi0<<std::endl;


}


///
/// @param n1 index for x1n1
/// @param t time
/// @return coefficient for evolution using H3
double os_DCE_Evolution::f(int n1, double t){

double x1n1Squared=x1ValsAllSquared[n1];

double val= -g0*omegac*std::sqrt(2.0/omegam)*std::sin(omegap*t)*x1n1Squared\
            +0.5*g0*std::sqrt(2.0/omegam)*std::sin(omegap*t);
    return val;


}


///
/// @param j time step
/// @param psi wavefunction at the beginning of the time step j
/// @return
arma::cx_dmat os_DCE_Evolution::evolution1Step(const int&j, const arma::cx_dmat& psi){
    arma::cx_dmat psiCurr(psi);

    double tj=timeValsAll[j];

    ///////////////////operator exp(-idt H1)
    //operator U15, for each column n2
    for (int n2=0;n2<N2;n2++){
        double x2n2=x2ValsAll[n2];
        psiCurr.col(n2)*=std::exp(1i*dt*0.5*g0*std::sqrt(2.0*omegam)*std::cos(omegap*tj)*x2n2);
    }





}