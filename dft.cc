#include <sndfile.h>

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

std::vector<double> ForReals(const std::vector<std::complex<double>>& v) {
  auto vw = std::views::transform(
      v, [](const std::complex<double>& c) { return c.real(); });
  return std::vector(vw.begin(), vw.end());
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

std::vector<double> GenerateRealSinusoid() {
  constexpr double kAmplitude = 0.8;
  constexpr double kFreqHz = 1000;
  double phi = std::numbers::pi / 2;

  auto time = arange(-0.002, 0.002, 1.0 / kSamplingRate);
  std::vector<double> signal;
  for (const auto& v : time) {
    signal.push_back(kAmplitude *
                     std::cos(2.0 * std::numbers::pi * kFreqHz * v + phi));
  }
  return signal;
}

std::vector<std::complex<double>> GenerateComplexSinusoid() {
  constexpr double k = 5;
  auto n = arange(-N / 2, N / 2);
  std::vector<std::complex<double>> signal;
  for (const auto& t : n) {
    double phase = 2.0 * std::numbers::pi * k * t / N;
    signal.push_back(std::polar(1.0, phase));
  }
  return signal;
}

std::vector<std::complex<double>> DFT(std::span<const double> signal) {
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

std::vector<std::complex<double>> iDFT(std::vector<std::complex<double>> dft) {
  std::vector<std::complex<double>> signal(dft.size());
  auto multp = 1.0 / dft.size();
  for (int n = 0; n < dft.size(); n++) {
    std::complex<double> accum(0, 0);
    for (int k = 0; k < dft.size(); k++) {
      double phase = 2.0 * std::numbers::pi * n * k / signal.size();
      std::complex<double> amp = std::polar(1.0, phase);
      accum += multp * dft[k] * amp;
    }
    signal[n] = accum;
  }
  return signal;
}

void OpenSquareWave() {
  constexpr std::string_view filename = "./square-wave-440.wav";

  SF_INFO fileinfo;
  SNDFILE* sndf = sf_open(filename.data(), SFM_READ, &fileinfo);

  std::cout << "YO - OPENED FILE:" << filename << " SR:" << fileinfo.samplerate
            << " Format:" << fileinfo.format
            << " Num frames:" << fileinfo.frames << std::endl;

  constexpr int NFFT = 512;
  constexpr int STEPSIZE = 256;
  const int num_frames = fileinfo.frames;
  const int num_channels = fileinfo.channels;

  std::vector<double> buffer(NFFT * num_channels, 0);
  int file_idx = 0;
  int write_num = 0;
  int num_read = 0;
  while (file_idx < num_frames) {
    int write_idx = 0;
    if (write_num == 1) {
      write_idx = buffer.size() / 2;
    }
    sf_count_t frames_read =
        sf_readf_double(sndf, &buffer[write_idx], STEPSIZE);
    std::cout << "REad: " << frames_read << std::endl;
    write_num++;
    if (write_num == 2) {
      // DO STFT
      buffer.clear();
    }
    file_idx += STEPSIZE;
  }
}

int main() {
  std::cout << "YO MO!" << std::endl;

  OpenSquareWave();

  // constexpr int kNumSamples = 1024;
  // auto signal =
  //     GenerateMultiToneSignal({{5440.0, 1.0}, {780.0, 0.5}, {1320.0, 0.25}},
  //                             kSamplingRate, kNumSamples);
  // auto time = arange(-0.002, 0.002, 1.0 / kSamplingRate);
  // auto time = arange(-N / 2, N / 2);
  // auto time = arange(0.0, 1.0, 1.0 / N);
  // //  auto signal = GenerateRealSinusoid();
  // auto signal = GenerateComplexSinusoid();
  // // GNUPlot(time, ForReals(signal));
  // GNUPlot(time, ForReals(signal));
  // auto X = DFT(ForReals(signal));
  // // auto X = DFT(signal);
  // auto [freqs, mags] = MagnitudeSpectrum(X, N);
  // GNUPlot(freqs, mags);

  // auto reconstructedSignal = iDFT(X);
  // GNUPlot(time, ForReals(reconstructedSignal));
}
