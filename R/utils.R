#' Reads out SHELL environment variable \code{MCCORES} and returns it as
#' integer. If the environment variable is not set, the value returned by
#' \code{detectCores()} is used as default. 
#'
#' @param env.var The name of the SHELL environment variable to read out.
#' Default is \code{MCCORES}.
#'
#' @export
#' @return An integer, either the value of the environment variable
#' \code{MCCORES} or the default value returned by
#' \code{parallel::detectCores()}.
getMcCores <- function(env.var = "MCCORES") {
    mc.cores <- as.integer(system("echo $MCCORES", intern = TRUE))
    if (!is.null(mc.cores) && !is.na(mc.cores) && mc.cores > 0) {
        message("Setting 'options(mc.cores = ", mc.cores, ")'")
        mc.cores
    } else {
        message("Setting 'options(mc.cores = ", detectCores(), ")'")
        detectCores()
    }
}
