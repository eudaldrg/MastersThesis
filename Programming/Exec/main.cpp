#include <iostream>

//#include <fftw3.h>
#include "fftw3.h"

#include "swift.h"

//C:\Program Files (x86)\mingw-w64\i686-8.1.0-posix-dwarf-rt_v6-rev0\mingw32\bin>g++ D:\GoogleDrive\TFM\Exec\main.cpp -o D:\GoogleDrive\TFM\Exec\a.exe -I "C:\Program Files (x86)\fftw" -I "C:\Program Files\boost\boost_1_71_0" -L "C:\Program Files (x86)\fftw" -l fftw3-3

void TestFFTW3()
{
    fftw_complex *in, *out;
    fftw_plan p;
    int N = 4;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    std::vector<std::complex<double>> x = {{1,0}, {2, -1}, {0, -1}, {-1, 2}};
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        in[i][0] = x[i].real();
        in[i][1] = x[i].imag();
    }
    fftw_execute(p); /* repeat as needed */
    std::vector<std::complex<double>> X = {{2,0}, {-2, -2}, {0, -2}, {4, 4}};
    std::cout << "Expected result";
    for (std::complex<double> val : X)
        std::cout << " " << val;
    std::cout << std::endl;
    std::cout << "Obtained result";
    for (std::size_t i = 0; i < N; ++i)
        std::cout << " (" << out[i][0] << "," << out[i][1];
    std::cout << std::endl;
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}

void BitReverseUT()
{
    {
        std::size_t x = 0b1101;
        std::size_t log2n = 4;
        std::size_t n = 0b1011;
        std::cout << "Reverse of x: " << Swift::bitReverse(x, log2n) << " expected " << n << std::endl;
    }
    {
        std::size_t x = 0b110101;
        std::size_t log2n = 6;
        std::size_t n = 0b101011;
        std::cout << "Reverse of x: " << Swift::bitReverse(x, log2n) << " expected " << n << std::endl;
    }
    {
        std::size_t x = 0b110101;
        std::size_t log2n = 4;
        std::size_t n = 0b1010;
        std::cout << "Reverse of x: " << Swift::bitReverse(x, log2n) << " expected " << n << std::endl;
    }
}

void FFTUT()
{
    std::vector<std::complex<double>> x = {{1,0}, {2, -1}, {0, -1}, {-1, 2}};
    std::vector<std::complex<double>> X = {{2,0}, {-2, -2}, {0, -2}, {4, 4}};
    std::vector<std::complex<double>> real_X = x;
    Swift::fft(real_X, false);
    std::cout << "Expected result";
    for (std::complex<double> val : X)
        std::cout << " " << val;
    std::cout << std::endl;
    std::cout << "Obtained result";
    for (std::complex<double> val : real_X)
        std::cout << " " << val;
    std::cout << std::endl;

    std::cout << "Original vec";
    for (std::complex<double> val : x)
        std::cout << " " << val;
    std::cout << std::endl;
    std::vector<std::complex<double>> real_x = real_X;
    Swift::fft(real_x, true);
    std::cout << "Recovered vec";
    for (std::complex<double> val : real_x)
        std::cout << " " << val;
    std::cout << std::endl;
}

int main() {
    TestFFTW3();
}
