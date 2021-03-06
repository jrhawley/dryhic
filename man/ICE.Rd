% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICE.R
\name{ICE}
\alias{ICE}
\title{Quick and dirty implementation of ICE matrix normalization}
\usage{
ICE(mat, itermax = 1, maxdev = 0.1, verbose = TRUE)
}
\arguments{
\item{mat}{A contact matrix (could be the output of \code{\link{get_contacts_matrix}})}

\item{maxdev}{Maximum relative deviation from the mean of the correction factor}

\item{verbose}{Logical indicating if progression messages should be turned on}

\item{intermax}{Maximum number of iterations}
}
\value{
A named \code{sparse_matrix} containing the normalized number of contacts per pair of genomic bins at the requested region plus the \code{b} attribute containing the vector of correction factors
}
\description{
This function takes a contact matrix (tipically the output of \code{\link{get_contacts_matrix}} and performs ICE normalization
}
\examples{
plot(0)
}
\references{
@article{imakaev2012iterative,
  title={Iterative correction of Hi-C data reveals hallmarks of chromosome organization},
  author={Imakaev, Maxim and Fudenberg, Geoffrey and McCord, Rachel Patton and Naumova, Natalia and Goloborodko, Anton and Lajoie, Bryan R and Dekker, Job and Mirny, Leonid A},
  journal={Nature methods},
  volume={9},
  number={10},
  pages={999--1003},
  year={2012},
  publisher={Nature Publishing Group}
}
}
