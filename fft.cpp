#include <iostream>
#include <complex>
#include <fftw3.h>
#include <vector>
#include <cmath>
#include <chrono>
using namespace std::complex_literals;
template<class T>
static void printVec(const std::vector <T> &vec) {
    for (int i = 0; i < vec.size() - 1; i++) {
        std::cout << vec[i] << ",";
    }
    std::cout << vec[vec.size() - 1] << std::endl;
}
int main() {

    int N1=3000;
    int N2=5000;

    std::complex<double>* x=new std::complex<double>[N1*N2];
    std::complex<double>* y=new std::complex<double>[N1*N2];

    for (int i=0;i<N1*N2;i++){
        double idb=static_cast<double >(i);
        x[i]=idb;
    }
    const auto tRowStart{std::chrono::steady_clock::now()};
    //row fft
    int rank=1;
    int n[]={N2};
    int howmany=N1;
    int istride=1;
    int ostride=1;
    int idist=N2;
    int odist=N2;
    int *inembed = n, *onembed = n;

    fftw_plan p= fftw_plan_many_dft(rank,n,howmany,reinterpret_cast<fftw_complex*>(x),inembed,
                                    istride,idist,reinterpret_cast<fftw_complex*>(y),
                                    onembed,ostride,odist,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    const auto tRowEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tRowEnd - tRowStart};
    std::cout<<"row fft time: "<< elapsed_secondsAll.count()  << " s" << std::endl;
    ////end of row fft

    //column fft
//    int rank=1;
//    int n[]={N1};
//    int howmany=N2;
//    int idist=1;
//    int odist=1;
//    int istride=N2;
//    int ostride=N2;
//    int *inembed = n, *onembed = n;
//
//    fftw_plan p= fftw_plan_many_dft(rank,n,howmany,reinterpret_cast<fftw_complex*>(x),inembed,
//                                    istride,idist,reinterpret_cast<fftw_complex*>(y),
//                                    onembed,ostride,odist,FFTW_FORWARD,FFTW_ESTIMATE);
//    fftw_execute(p);
//    fftw_destroy_plan(p);



    //////end of column fft

//    for(int i=0;i<N1;i++){
//        for(int j=0;j<N2;j++){
//            std::cout<<y[i*N2+j]<<",";
//        }
//        std::cout<<std::endl;
//    }

    delete []x;
    delete []y;

    return 0;
}
