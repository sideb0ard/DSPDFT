#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <numbers>
#include <ranges>
#include <vector>

constexpr double N = 500;
constexpr double kSamplingRate = 44100;

template <typename T>
std::vector<T> arange(T start = 0, T stop = N, T step = 1) {
  std::vector<T> values;
  for (T value = start; value < stop; value += step) {
    values.push_back(value);
  }
  return values;
}

struct Tone {
  double freq_hz;
  double amplitude;
};

std::vector<double> GenerateMultiToneSignal(const std::vector<Tone>& tones,
                                            double sampling_rate,
                                            int num_samples) {
  std::vector<double> signal(num_samples, 0.0);
  for (int n = 0; n < num_samples; n++) {
    double t = n / sampling_rate;
    double sample = 0.0;
    for (const auto& tone : tones) {
      sample +=
          tone.amplitude * std::cos(2.0 * std::numbers::pi * tone.freq_hz * t);
    }
    signal[n] = sample;
  }
  return signal;
}

std::pair<std::vector<double>, std::vector<double>> MagnitudeSpectrum(
    const std::vector<std::complex<double>>& X, double sampling_rate) {
  const int N = static_cast<int>(X.size());
  const int half = N / 2;

  std::vector<double> freqs(half);
  std::vector<double> mags(half);
  for (int k = 0; k < half; k++) {
    freqs[k] = k * sampling_rate / N;
    mags[k] = std::abs(X[k]) / N;
  }
  return {freqs, mags};
}

auto realView(const std::vector<std::complex<double>>& v) {
  return std::views::transform(
      v, [](const std::complex<double>& c) { return c.real(); });
}

void GNUPlot(const std::vector<double>& time, const std::vector<double>& amp) {
  struct PipeCloser {
    void operator()(FILE* f) const {
      if (f) pclose(f);
    }
  };
  using PipePtr = std::unique_ptr<FILE, PipeCloser>;
  PipePtr gp(popen("gnuplot -persist", "w"));
  if (!gp) {
    fprintf(stderr, "Could not open pipe to gnuplot\n");
    return;
  }
  fprintf(gp.get(), "set title 'Signal'\n");
  fprintf(gp.get(), "set xlabel 'Time (s)'\n");
  fprintf(gp.get(), "set ylabel 'Amplitude'\n");
  fprintf(gp.get(), "plot '-' with lines title 'signal'\n");

  for (size_t i = 0; i < time.size(); i++) {
    fprintf(gp.get(), "%f %f\n", time[i], amp[i]);
  }
  fprintf(gp.get(), "e\n");  // end of inlien data

  // if (pclose(gp) == -1) {
  //   fprintf(stderr, "uh oh\n");
  // }
}

void GenerateRealSinusoid() {
  constexpr double kAmplitude = 0.8;
  constexpr double kFreqHz = 1000;
  double phi = std::numbers::pi / 2;

  auto time = arange(-0.002, .002, 1.0 / kSamplingRate);
  std::cout << "GOT VALS: len :" << time.size() << "\n";
  for (const auto& v : time) {
    std::cout << v << std::endl;
  }

  std::vector<double> signal;
  for (const auto& v : time) {
    signal.push_back(kAmplitude *
                     std::cos(2.0 * std::numbers::pi * kFreqHz * v + phi));
  }

  GNUPlot(time, signal);
}

void GenerateComplexSinusoid() {
  constexpr double k = 3;
  auto n = arange(-N / 2, N / 2);
  std::vector<std::complex<double>> signal;
  for (const auto& t : n) {
    double phase = 2.0 * std::numbers::pi * k * t / N;
    signal.push_back(std::polar(1.0, phase));
  }
  auto view = realView(signal);
  GNUPlot(n, std::vector<double>(view.begin(), view.end()));
}

std::vector<std::complex<double>> DFT(const std::vector<double>& signal) {
  std::vector<std::complex<double>> X(signal.size());
  for (int k = 0; k < signal.size(); k++) {
    std::complex<double> accum(0, 0);
    for (int n = 0; n < signal.size(); n++) {
      double phase = 2.0 * std::numbers::pi * n * k / signal.size();
      std::complex<double> spectra = std::polar(1.0, phase);
      accum += signal[n] * std::conj(spectra);
    }
    X[k] = accum;
  }
  return X;
}

int main() {
  std::cout << "YO MO!" << std::endl;
  constexpr int kNumSamples = 1024;
  auto signal =
      GenerateMultiToneSignal({{5440.0, 1.0}, {780.0, 0.5}, {1320.0, 0.25}},
                              kSamplingRate, kNumSamples);
  auto X = DFT(signal);
  auto [freqs, mags] = MagnitudeSpectrum(X, kSamplingRate);
  GNUPlot(freqs, mags);
}
