\name{noc_fit}
\alias{noc_fit}
\title{Fit a member of the contextual empirical Bayes curve family}
\usage{noc_fit(y, method = c("EB0", "RAW", "NOC"))}
\arguments{
  \item{y}{Finite numeric observations in equally spaced design order, with
  at least eight points.}
  \item{method}{One or more family members. The default fits all three
  with the same fitted levelwise base priors.}
}
\value{
A list with the prediction matrix \code{prediction}, with one column per
requested member, the named vector \code{strength} of selected strengths
in \eqn{[0,1]}, and the difference-based scale estimate \code{sigma}.
}
\description{
Uses a fixed stationary wavelet transform and shared levelwise normal
mixture priors. EB0 fixes context strength at zero. RAW and NOC select it
with a penalized composite marginal criterion. NOC projects the context
directions off the target analysis row before measuring context energy.
}
\details{
The procedure assumes independent, homoscedastic Gaussian noise on the original
equally spaced observations. For non-dyadic lengths the transform uses
the fixed source-copy adapter. Copied observations do not increase the
sample size. A constant input is returned exactly. Nonconstant inputs with
zero estimated scale or unsupported numerical geometry raise errors.
}
