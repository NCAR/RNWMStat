#Calculate the 'annual mean of the flow exceeded 10% of the time' (Boscarello paper)

#-------FUNCTION SUMMARY -----------
#1. Q10_50
#   REQUIRED INPUTS
#      1. Spline function fitted to the flow
#   OPTIONAL INPUTS
#       None
#   DEPENDENCIES
#       None
#-----------------------------------

Q10_50 <- function (splineFlow, na.rm=TRUE) {
    Q10.50 <- splineFlow(0.1)/splineFlow(0.5)
    return(Q10.50)
}
