// Bit-identical C++ port of ASCAT's aspcfpart() (the O(N^2) core of fastAspcf).
// The accumulation order, the round/which.min semantics and the trace-back are
// replicated exactly, so results match the pure-R version to the last bit.
#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector aspcfpart_cpp(NumericVector logRpart, NumericVector allBflip,
                            double a, double b, double d,
                            double sd1, double sd2, double Nfull,
                            int kmin, double gamma) {
  double from    = std::max(1.0, a);
  double usefrom = std::max(1.0, a + d);
  double useto   = std::min(Nfull, b - d);

  int N = logRpart.size();
  if (N < 2 * kmin) return NumericVector::create(0.0);

  // R is 1-based and the algorithm is written that way; keep the same indexing
  // by copying the inputs into 1-based buffers (element 0 unused).
  std::vector<double> y1(N + 1, 0.0), y2(N + 1, 0.0);
  for (int i = 1; i <= N; ++i) { y1[i] = logRpart[i - 1]; y2[i] = allBflip[i - 1]; }

  // 1-based scratch arrays (index 0 unused) to mirror the R code exactly
  std::vector<double> bestCost(N + 1, 0.0), bestSplit(N + 1, 0.0);
  std::vector<double> bestAver1(N + 1, 0.0), bestAver2(N + 1, 0.0);
  std::vector<double> Sum1(N + 1, 0.0), Sum2(N + 1, 0.0);
  std::vector<double> Kvad1(N + 1, 0.0), Kvad2(N + 1, 0.0);
  std::vector<double> Aver1(N + 1, 0.0), Aver2(N + 1, 0.0);
  std::vector<double> Cost(N + 1, 0.0);

  double initSum1 = 0.0, initKvad1 = 0.0, initSum2 = 0.0, initKvad2 = 0.0;
  for (int i = 1; i <= kmin; ++i) {
    initSum1  += y1[i]; initKvad1 += y1[i] * y1[i];
    initSum2  += y2[i]; initKvad2 += y2[i] * y2[i];
  }
  double initAve1 = initSum1 / kmin;
  double initAve2 = initSum2 / kmin;
  double sd1sq = sd1 * sd1, sd2sq = sd2 * sd2;

  bestCost[kmin]  = (initKvad1 - initSum1 * initAve1) / sd1sq
                  + (initKvad2 - initSum2 * initAve2) / sd2sq;
  bestAver1[kmin] = initAve1;
  bestAver2[kmin] = initAve2;

  int kminP1 = kmin + 1;
  for (int k = kminP1; k <= 2 * kmin - 1; ++k) {
    double yk1 = y1[k], yk2 = y2[k];
    for (int j = kminP1; j <= k; ++j) {
      Sum1[j]  += yk1;
      Aver1[j]  = Sum1[j] / (double)(k - j + 1);
      Kvad1[j] += yk1 * yk1;
      Sum2[j]  += yk2;
      Aver2[j]  = Sum2[j] / (double)(k - j + 1);
      Kvad2[j] += yk2 * yk2;
    }
    bestAver1[k] = (initSum1 + Sum1[kminP1]) / k;
    bestAver2[k] = (initSum2 + Sum2[kminP1]) / k;
    // R parses `k*bestAver1[k]^2` as k * (x*x) -- keep that association exactly
    double cost1 = ((initKvad1 + Kvad1[kminP1]) - k * (bestAver1[k] * bestAver1[k])) / sd1sq;
    double cost2 = ((initKvad2 + Kvad2[kminP1]) - k * (bestAver2[k] * bestAver2[k])) / sd2sq;
    bestCost[k] = cost1 + cost2;
  }

  for (int n = 2 * kmin; n <= N; ++n) {
    int nMkminP1 = n - kmin + 1;
    double yn1 = y1[n], yn2 = y2[n];
    for (int j = kminP1; j <= n; ++j) {
      Sum1[j]  += yn1;
      Aver1[j]  = Sum1[j] / (double)(n - j + 1);
      Kvad1[j] += yn1 * yn1;
      Sum2[j]  += yn2;
      Aver2[j]  = Sum2[j] / (double)(n - j + 1);
      Kvad2[j] += yn2 * yn2;
    }
    double cost1 = 0.0, cost2 = 0.0;
    for (int j = kminP1; j <= nMkminP1; ++j) {
      cost1 = (Kvad1[j] - Sum1[j] * Aver1[j]) / sd1sq;
      cost2 = (Kvad2[j] - Sum2[j] * Aver2[j]) / sd2sq;
      Cost[j] = bestCost[j - 1] + cost1 + cost2;   // bestCost[kmin:(n-kmin)], kmin = j-1
    }
    // which.min(): first index attaining the minimum
    int pos = kminP1;
    for (int j = kminP1 + 1; j <= nMkminP1; ++j) if (Cost[j] < Cost[pos]) pos = j;
    double cost  = Cost[pos] + gamma;
    double aver1 = Aver1[pos], aver2 = Aver2[pos];

    double totAver1 = (Sum1[kminP1] + initSum1) / n;
    double totCost1 = ((Kvad1[kminP1] + initKvad1) - n * totAver1 * totAver1) / sd1sq;
    double totAver2 = (Sum2[kminP1] + initSum2) / n;
    double totCost2 = ((Kvad2[kminP1] + initKvad2) - n * totAver2 * totAver2) / sd2sq;
    double totCost  = totCost1 + totCost2;

    if (totCost < cost) { pos = 1; cost = totCost; aver1 = totAver1; aver2 = totAver2; }
    bestCost[n] = cost; bestAver1[n] = aver1; bestAver2[n] = aver2;
    bestSplit[n] = (double)(pos - 1);
  }

  // trace back, prepending exactly like the R while-loop
  std::vector<double> bp;
  bp.push_back((double)N);
  double nn = (double)N;
  while (nn > 0) {
    double sp = bestSplit[(int)nn];
    bp.push_back(sp);
    nn = sp;
  }
  std::reverse(bp.begin(), bp.end());

  std::vector<double> out;
  out.reserve(bp.size());
  for (size_t i = 0; i < bp.size(); ++i) {
    double v = bp[i] + from - 1;
    if (v >= usefrom && v <= useto) out.push_back(v);
  }
  return wrap(out);
}
