#include <iostream>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <matplot/matplot.h>
#include "AudioFile.h"
#include <complex>
#include <cmath>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

void plot(py::array_t<double> t, py::array_t<double> x){
    std::vector<double> time(t.size());
    std::vector<double> signal(x.size());

    auto buf_t = t.request();
    auto buf_x = x.request();

    double* ptr_t = static_cast<double*>(buf_t.ptr);
    double* ptr_x = static_cast<double*>(buf_x.ptr);

    std::memcpy(time.data(), ptr_t, buf_t.size * sizeof(double));
    std::memcpy(signal.data(), ptr_x, buf_x.size * sizeof(double));

    auto fig = matplot::figure();
    matplot::title("Podany sygnal:");
    matplot::plot(time, signal);
    matplot::show();
}

void grafuj(const std::string& filePath) {
     AudioFile<double> audioFile;
   audioFile.load(filePath);

    std::vector<std::vector<double>> samples = audioFile.samples;
    std::vector<double> mono = samples[0];
     if (mono.size() > 1000) {
        mono.resize(1000);
    }
    auto fig = matplot::figure();
    matplot::title("Podany sygnal:");
    matplot::plot(mono);
    matplot::show();
}

std::vector<std::complex<double>> dft(const std::vector<std::complex<double>>& samples) {
    const size_t N = samples.size();
    std::vector<std::complex<double>> result(N);
    std::complex<double> M_I_2PI_DL = (std::complex<double>(0, -2.0 * M_PI) / static_cast<double>(samples.size()));

    for (size_t i = 0; i < N; i++) {
        result[i] = std::complex<double>(0, 0);
        for (size_t n = 0; n < N; n++) {
            result[i] += samples[n] * std::exp(M_I_2PI_DL * (double)i * (double)n);
        }
    }

    return result;
}

void progowanie(double prog, const std::string& filePath){
    AudioFile<double> audioFile;
    audioFile.load(filePath);

    std::vector<std::vector<double>> samples = audioFile.samples;
    std::vector<double> mono = samples[0];
    if (mono.size() > 1000) {
        mono.resize(1000);
    }

    int N = mono.size();

    std::vector<int> progowany(N, 0);
    for(int i = 0; i<N; i++){
        if(mono[i] > prog) progowany[i] = 1;
        else progowany[i] = 0;
    }

    auto fig = matplot::figure();
    matplot::plot(progowany);
    matplot::title("Sygnal progowany:");
    matplot::show();
}

void transformata(const std::string& filePath) {
 AudioFile<double> audioFile;
    audioFile.load(filePath);

    std::vector<std::vector<double>> samples = audioFile.samples;
    std::vector<double> mono = samples[0];
    if (mono.size() > 1000) {
        mono.resize(1000);
    }

    int N = mono.size();


    std::vector<std::complex<double>> complex_samples(N);
    for (int i = 0; i < N; i++) {
        complex_samples[i] = std::complex<double>(mono[i], 0.0);
    }


    std::vector<std::complex<double>> dft_result = dft(complex_samples);

    std::vector<double> amplitudes(N);
    for (int i = 0; i < N; ++i) {
        amplitudes[i] = std::abs(dft_result[i]);
    }

    auto fig = matplot::figure();
    matplot::plot(amplitudes);
    matplot::title("Dyskretna Transformata Fouriera podanego sygnału:");
    matplot::show();
}

std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& samples) {
    const size_t N = samples.size();
    std::vector<std::complex<double>> result(N);
    std::complex<double> M_I_2PI_DL = (std::complex<double>(0, 2.0 * M_PI) / static_cast<double>(samples.size()));

    for (size_t i = 0; i < N; i++) {
        result[i] = std::complex<double>(0, 0);
        for (size_t n = 0; n < N; n++) {
            result[i] += samples[n] * std::exp(M_I_2PI_DL * (double)i * (double)n);
        }
            result[i] /= static_cast<double>(samples.size());
    }

    return result;
}

void odwrotna(const std::string& filePath){
AudioFile<double> audioFile;
    audioFile.load(filePath);

    std::vector<std::vector<double>> samples = audioFile.samples;
    std::vector<double> mono = samples[0];
    if (mono.size() > 1000) {
        mono.resize(1000);
    }

    int N = mono.size();


    std::vector<std::complex<double>> complex_samples(N);
    for (int i = 0; i < N; i++) {
        complex_samples[i] = std::complex<double>(mono[i], 0.0);
    }


    std::vector<std::complex<double>> idft_result = idft(complex_samples);

    std::vector<double> invAmplitudes(N);
    for (int i = 0; i < N; ++i) {
        invAmplitudes[i] = std::abs(idft_result[i]);
    }

    auto fig = matplot::figure();
    matplot::plot(invAmplitudes);
    matplot::title("Odwrotna Dyskretna Transformata Fouriera podanego sygnału:");
    matplot::show();
}
namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";
    m.def("transformata", &transformata);
    m.def("grafuj", &grafuj);
    m.def("dft", &dft);
    m.def("idft", &idft);
    m.def("odwrotna", &odwrotna);
    m.def("plot", &plot);
    m.def("progowanie", &progowanie);
    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
