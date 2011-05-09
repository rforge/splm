\name{spml}
\alias{spml}
\title{Spatial Panel Model by Maximum Likelihood}
\description{Maximum likelihood (ML) estimation of spatial panel models, possibly with fixed or random effects.}


\usage{
spml(formula, data, index = NULL, listw, listw2=listw,
                 model=c("within","random","pooling"),
                 effect=c("individual","time","twoways"),
                 lag=FALSE, spatial.error=c("b","kkp","none"),
                 ...)
}

\arguments{
  \item{formula}{a symbolic description of the model to be estimated}
  \item{data}{an object of class \code{\link{data.frame}} or \code{pdata.frame}. A data frame containing the variables
  in the model. When the object is a \code{\link{data.frame}}, the first two columns shall contain the indexes, unless otherwise specified. See \code{index}}
  \item{index}{if not NULL (default), a character vector to identify the indexes among the columns of the \code{\link{data.frame}}}
  \item{listw}{an object of class \code{listw} or a \code{matrix}. It represents the spatial weights to be used in estimation.}
  \item{listw2}{an object of class \code{listw} or a
  \code{matrix}. Second of set spatial weights for estimation, if
  different from the first (e.g., in a 'sarar' model).}
  \item{model}{one of \code{c("within", "random", "pooling").}}
  \item{effect}{one of \code{c("individual","time","twoways")}; the
  effects introduced in the model.}
  \item{lag}{default=\code{FALSE}. If \code{TRUE}, a spatial lag of the dependent variable is added.}
  \item{spatial.error}{one of \code{c("b","kkp","none")}. The type of
  spatial error in the specification, if any. See details.}
  \item{...}{additional argument to pass over to other functions}
}
\details{
  The models are estimated by two-step Maximum Likelihood.
  Further optional parameters to be passed on to the estimator may be:
pvar: if \code{TRUE} the \code{pvar} function is called
hess: if \code{TRUE} use numerical Hessian instead of GLS for the
standard errors of the estimates
quiet: if \code{FALSE} report function and parameters values during
optimization
initval: one of \code{c("zeros", "estimate")}, the initial values for
  the parameters. If \code{"zeros"} a vector of zeros is used. if
  \code{"estimate"} the initial values are retreived from the estimation
  of the nested specifications. Alternatively, a numeric vector can be
  specified.
  x.tol: Tolerance. See \code{\link{nlminb}} for details.
  rel.tol: Relative tolerance. See \code{\link{nlminb}} for details.
}
\value{
  An object of class \code{"splm"}.
  \item{coefficients}{coefficients estimate of the model parameters }
  \item{arcoef}{the coefficient for the spatial lag on \code{y}}
  \item{errcomp}{the estimates of the error variance components}
  \item{vcov}{the asymptotic variance covariance matrix of the estimated coefficients}
  \item{vcov.arcoef}{the asymptotic variance of the
  estimated spatial lag parameter}
  \item{vcov.errcomp}{the asymptotic variance covariance matrix of the
  estimated error covariance parameters}
  \item{type}{'random effects ML'}
  \item{residuals}{the model residuals}
  \item{fitted.values}{the fitted values, calculated as \eqn{\hat{y}=X \hat{\beta}}}
  \item{sigma2}{GLS residuals variance}
  \item{model}{the matrix of the data used}
  \item{call}{the call used to create the object}
  \item{logLik}{the value of the log likelihood function at the optimum }
  \item{errors}{the value of the \code{errors} argument}
}
\references{
Baltagi, B.H., Song, S.H., Jung B. and Koh, W. (2007)
Testing panel data regression models with spatial and serial error correlation.
\emph{Journal of Econometrics}, \bold{140}, 5-51.
}

\author{Giovanni Millo}

\seealso{\code{\link{spgm}}, \code{\link{spfeml}}}

\examples{
data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1974, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
## the two standard specifications:
## random effects panel with spatial errors
respaterr <- spml(fm, data = Produc, listw = mat2listw(usaww), model="random", spatial.error="b")
summary(respaterr)
## random effects panel with spatial lag
respatlag <- spml(fm, data = Produc, listw = mat2listw(usaww), model="random", spatial.error="b", lag=TRUE)
summary(respatlag)
}

\keyword{spatial}