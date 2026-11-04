import math
import statistics

def mean(values):
    return math.fsum(values) / len(values)


def mcse(values):
    return math.sqrt(statistics.variance(values) / len(values)) if len(values) > 1 else None


def losses(errors):
    assert errors and all(math.isfinite(x) for x in errors)
    return dict(MSE=mean([x*x for x in errors]), MAE=mean([abs(x) for x in errors]),
                MAXAE=max(abs(x) for x in errors))


def decomposition(errors):
    R, n = len(errors), len(errors[0])
    assert R and n and all(len(e) == n for e in errors)
    point_means = [mean([e[i] for e in errors]) for i in range(n)]
    apparent = mean([v*v for v in point_means])
    if R == 1:
        return dict(observed_squared_mean_error=apparent, corrected_squared_bias=None,
                    mean_sample_variance=None, identity_residual=None)
    variance = mean([statistics.variance([e[i] for e in errors]) for i in range(n)])
    corrected = apparent - variance/R
    residual = mean([losses(e)['MSE'] for e in errors]) - corrected - variance
    assert abs(residual) <= 1e-12
    return dict(observed_squared_mean_error=apparent, corrected_squared_bias=corrected,
                mean_sample_variance=variance, identity_residual=residual)


def stratified_mcse(strata):
    if any(len(v) < 2 for v in strata):
        return None
    return math.sqrt(math.fsum(statistics.variance(v)/len(v) for v in strata))/len(strata)


