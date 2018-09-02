#' Hourly measurements of 6 air pollutant concentrations in Bloomsbury, London,
#' UK.
#'
#' A dataset containing the concentration of 6 pollutants in Bloomsbury, London,
#' UK with hourly measurements from 1st January 2000 to 3st December 2017.
#' Courtesy of King's College Air Pollution monitoring network. Values expressed
#' in micrograms per cubic meter.
#'
#' @docType data
#' @usage data(hourly_bloomsbury_air_pollution_2000_2017)
#' @format A data frame with 157710 rows and 9 variables: \describe{
#'   \item{index}{index of entry}{date of entry in format
#'   dd/mm/yyyy}
#'   \item{time}{time of day of entry in format HH:MM}
#'   \item{O3}{concentration in Ozone, in mu g / m^3}
#'   \item{NO}{concentration in
#'   Nitrogen Oxide, in mu g / m^3}
#'   \item{NO2}{concentration in Nitrogen
#'   Dioxyde, in mu g / m^3} \item{PM10}{concentration in particulate matter 10
#'   micrometers or less in diameter, in mu g / m^3}
#'   \item{SO2}{concentration in
#'   Sulfur Dioxide, in mu g / m^3}}
#' @source \url{https://www.londonair.org.uk/} data selection tool.
#' @keywords datasets
#' @examples
#' data(hourly_bloomsbury_air_pollution_2000_2017)
#' \donttest{plot(hourly_bloomsbury_air_pollution_2000_2017$O3)}
"hourly_bloomsbury_air_pollution_2000_2017"
