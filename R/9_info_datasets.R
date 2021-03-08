#############################################
#' Attendance Records data set
#' @description School administrators study the attendance behavior of high school juniors at two schools. 
#' @docType data
#' @usage data("Attendance")
#' @format Data frame containing 314 observations on 4 variables.
#' \describe{
#'   \item{daysabs}{number of days absent.}
#'   \item{gender}{gender of the student.}
#'   \item{prog}{three-level factor indicating the type of instructional program in which the student is enrolled.}
#'   \item{math}{standardized math score.}
#' }
#' 
#' @details 
#' School administrators study the attendance behavior of high school juniors at two schools. Predictors of the number of days of absence include the type of program in which 
#' the student is enrolled and a standardized test in math. Attendance data on 314 high school juniors from two urban high schools. 
#' The response variable of interest is days absent, \code{daysabs}. The variable \code{math} is the standardized math score for each student. 
#' The variable \code{prog} is a three-level factor indicating the type of instructional program in which the student is enrolled.
#' 
#' @references
#' Hughes, M. and Fisher, T. (2020)  \href{http://www.users.miamioh.edu/fishert4/sta363/}{Introduction to Statistical Modeling.}
#'
#' @source Data can be obtained from \href{https://github.com/tjfisher19/introStatModeling}{Introduction to Statistical Modeling Github Repository}. See also \emph{Barreto-Souza and Simas (2020)} for further details.
#' @examples
#' data("Attendance", package = "mixpoissonreg")
#' 
#' daysabs_fit <- mixpoissonregML(daysabs ~ gender + math + prog | gender + 
#' math + prog, data = Attendance)
#' summary(daysabs_fit)
#' 
"Attendance"