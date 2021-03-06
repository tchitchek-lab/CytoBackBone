% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{FCS-class}
\alias{FCS-class}
\alias{FCS}
\title{FCS class definition}
\description{
FCS are S4 objects containing marker cell expressions obtained from cytometry profiles.
}
\details{
This object mainly stores for each cell the intensities of all cell markers.

The slot 'trans.para' is a named list contains different parameters depending of the transformation applied on the marker expression intensities. The scale (cofactor) of the arcsinh transformation function is parametrized using the 'arcsinh.scale' value. The shift of the log transformation function is parametrized using the 'log.shift' value and the base of the log transformation function is parametrized using the 'log.base' value. If no transformation function have been applied, the 'trans.para' slot is set to NULL.
}
\section{Slots}{

\describe{
\item{\code{name}}{a character indicating the internal name of the FCS object}

\item{\code{profiles}}{a character vector containing the names of the FCS profiles}

\item{\code{cell.nb}}{an integer value indicating the number of FCS profiles}

\item{\code{markers}}{a character vector containing the marker names}

\item{\code{markers.nb}}{an integer value indicating the number of markers}

\item{\code{intensities}}{a numeric matrix containing the intensities of each marker for each FCS profile}

\item{\code{trans}}{a character specifying the name of a transformation function applied on the marker expression intensities. Possible values are "arcsinh" for arc sin hyperbolic transformation, "log" for logarithmic transformation, or "none" for no transformation}

\item{\code{trans.para}}{a named list containing parameters of the transformation. Please refer to the details section for more details}

\item{\code{trans.exclude}}{a character vector containing the marker names for which no transformation has been applied on}
}}

