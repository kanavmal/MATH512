# Evaluate the following expected value and probability:
# E[X_2^0.6] and P(X_2 > 2)
# Where the Ito's Processes evolve according to the following stochastic differential equation:
# dX_t = (1/4 + 1/3X_t)dt + 3/5dW_t, X_0 = 2
# and W_t is a standard Wiener process.

using Random
using Statistics
using Plots
using LaTeXStrings
using DifferentialEquations
 



# Save the plot
savefig("imgs/question3_plot.png")