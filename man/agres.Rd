\name{agres}
\alias{agres}
\docType{data}
\title{European agricultural data}
\description{Data on research expenditure, research activity, productivity and impacts of Agriculture from 1980 to 2014 in EU 15 countries. Source: EUROSTAT, Faostat, OECD, USDA.}
\usage{data(agres)}
\format{
  A data frame with a total of 490 observations on the following 19 variables:
  \describe{
    \item{\code{COUNTRY}}{Code of the country. Belgium and Luxembourg are merged together into the code \code{BE}.}
    \item{\code{YEAR}}{Year.}
    \item{\code{GDP}}{Gross domestic product (international dollars).}
    \item{\code{EMPL_AGR}}{Persons employed in Agriculture (count).}
    \item{\code{UAA}}{Utilized agricultural area (hectares).}
    \item{\code{PATENT_OTHER}}{Mechanical, chemical and environment-related patent applications (count).}
    \item{\code{GBAORD_AGR}}{Government research expenditure (million euro PPS).}
    \item{\code{BERD_AGR}}{Business enterprise research expenditure (million euro PPS).}
    \item{\code{RD_EDU_AGR}}{Agricultural researchers with tertiary education (count).}
    \item{\code{EU_PRO_AGR}}{Projects funded by the European Union (count).}
    \item{\code{PATENT_AGR}}{Agricultural patent applications (count).}
    \item{\code{TFPC}}{Total factor productivity (index 2005=100).}
    \item{\code{EPC}}{Environmental performance (index 2005=100).}
    \item{\code{GVA_AGR}}{Gross value added of Agriculture (international dollars).}
    \item{\code{PPI_AGR}}{Producer price of agricultural output (index 2004=100).}
    \item{\code{ENERGY_RENEW}}{Share of renewable energy in gross final energy consumption (index 2004=100).}
    \item{\code{INCOME_RURAL}}{Mean familiar income in rural areas (international dollars).}
    \item{\code{UNEMPL_RURAL}}{Unemployment rate in rural areas (index 2004=100).}
    \item{\code{HEALTH_RURAL}}{Quota of persons with at least a good health status in rural areas (index 2004=100).}
  }
}
