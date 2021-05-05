#include <iostream>
#include <complex>
#include <cmath>
#include <vector>

std::vector<std::complex<double> > dft(std::vector<std::complex<double> > data){
    /*
     return dft of data;
     data should contain N samples which were taken from one period
     */
    int N = data.size();
    std::vector<std::complex<double> > output;
    output.reserve(N);
    for(int k=0; k < N; ++k){
        std::complex<double> Sum = std::complex<double>(0.0, 0.0);
        for(int n=0; n < N; ++n){
            double realPart = cos(2 * M_PI / N * k * n);
            double imPart = -sin(2 * M_PI / N * k * n);
            Sum += data[n] * std::complex<double>(realPart, imPart);
        }
        Sum /= std::complex<double> (N, 0);
        output.push_back(Sum);
    }
    return output;
}

int main(){
    int N = 1024;
    std::vector<std::complex<double> > data(N);
    int f = 3;
    for(int i = 0; i < N; i++){
        data[i] = std::complex<double>(cos(2 * M_PI * f * i / N), 0.0);
    }
    std::vector<std::complex<double> > ans = dft(data);
    for(int i = 1; i < 5; i++){
        std::cout << ans[i].real() << " " << ans[i].imag() << '\n';
    }
    return 0;
}
//
// Created by Arseniy on 05.05.2021.
//

