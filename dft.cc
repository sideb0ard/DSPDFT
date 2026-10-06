#include <sndfile.h>

#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <numbers>
#include <ranges>
#include <vector>

constexpr double kComplexSinusoidLength = 500;
constexpr double kSamplingRate = 44100;
constexpr int kNfft = 512;
constexpr int kStep = 256;

template <typename T>
std::vector<T> arange(T start = 0, T stop = kComplexSinusoidLength,
                      T step = 1) {
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
void GNUSPlot(const std::vector<std::vector<double>>& heatmap,
              double time_step_sec, double bin_hz) {
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
  fprintf(gp.get(), "set pm3d map\n");
  fprintf(gp.get(), "set palette rgbformulae 22,13,-31\n");
  fprintf(gp.get(), "set title 'Heatmap'\n");
  fprintf(gp.get(), "set xlabel 'Time (s)'\n");
  fprintf(gp.get(), "set ylabel 'Frequency (Hz)'\n");
  fprintf(gp.get(), "splot '-' using 1:2:3 with pm3d notitle\n");

  for (size_t m = 0; m < heatmap.size(); ++m) {
    const double t = m * time_step_sec;
    for (size_t k = 0; k < heatmap[m].size(); ++k) {
      const double f = k * bin_hz;
      fprintf(gp.get(), "%f %f %f\n", t, f, heatmap[m][k]);
    }
    fprintf(gp.get(), "\n");
  }
  fprintf(gp.get(), "e\n");
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
  auto n = arange(-kComplexSinusoidLength / 2, kComplexSinusoidLength / 2);
  std::vector<std::complex<double>> signal;
  for (const auto& t : n) {
    double phase = 2.0 * std::numbers::pi * k * t / kComplexSinusoidLength;
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

std::vector<double> OpenWaveFile(const std::string_view filename) {
  SF_INFO fileinfo;
  SNDFILE* sndf = sf_open(filename.data(), SFM_READ, &fileinfo);

  const int num_frames = fileinfo.frames;
  const int num_channels = fileinfo.channels;

  std::cout << "YO - OPENED FILE:" << filename << " SR:" << fileinfo.samplerate
            << " Format:" << fileinfo.format << " Num frames:" << num_frames
            << std::endl;

  std::vector<double> buffer(num_frames * num_channels, 0);
  sf_count_t frames_read = sf_readf_double(sndf, buffer.data(), num_frames);
  sf_close(sndf);

  return buffer;
}

void FFT(std::vector<std::complex<double>>& x) {
  const size_t N = x.size();
  if (N <= 1) return;

  std::vector<std::complex<double>> even(N / 2);
  std::vector<std::complex<double>> odd(N / 2);
  for (size_t i = 0; i < N / 2; ++i) {
    even[i] = x[2 * i];
    odd[i] = x[2 * i + 1];
  }
  FFT(even);
  FFT(odd);

  for (size_t k = 0; k < N / 2; k++) {
    std::complex<double> twiddle =
        std::polar(1.0, -2.0 * M_PI * k / N) * odd[k];
    x[k] = even[k] + twiddle;
    x[k + N / 2] = even[k] - twiddle;
  }
}

void iFFT(std::vector<std::complex<double>>& x) {
  const size_t n = x.size();

  // conjugate input
  for (auto& v : x) v = std::conj(v);

  FFT(x);

  // conjugate result and normalize ny N
  for (auto& v : x) v = std::conj(v) / static_cast<double>(n);
}

void STFT() {
  auto chirp = OpenWaveFile("chirp.wav");

  std::vector<double> window(kNfft, 0);
  for (int n = 0; n < kNfft; n++) {
    // double hann = 0.5 - 0.5 * std::cos(2.0 * M_PI * n / kNfft);
    // window[n] = std::sqrt(hann);
    window[n] = 0.5 - 0.5 * std::cos(2.0 * M_PI * n / kNfft);
  }

  const size_t num_frames = chirp.size() >= static_cast<size_t>(kNfft)
                                ? (chirp.size() - kNfft) / kStep + 1
                                : 0;

  std::vector<std::vector<double>> spectrogram(
      num_frames, std::vector<double>(kNfft / 2 + 1));

  for (size_t m = 0; m < num_frames; m++) {
    const size_t start = m * kStep;
    std::vector<std::complex<double>> frame(kNfft);
    for (int n = 0; n < kNfft; n++) {
      frame[n] = std::complex<double>(chirp[start + n] * window[n], 0.0);
    }

    FFT(frame);

    for (int k = 0; k < kNfft / 2 + 1; ++k) {
      double mag = std::abs(frame[k]);
      spectrogram[m][k] = 20.0 * std::log10(mag + 1e-9);
    }
  }
  constexpr double kFrameSec = kStep / kSamplingRate;
  constexpr double kBinHz = kSamplingRate / kNfft;
  GNUSPlot(spectrogram, kFrameSec, kBinHz);
}

int main() {
  std::cout << "YO MO!" << std::endl;

  STFT();

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
