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
 
# Define the SDE
function f(du, u, p, t)
    du[1] = (1/4 + 1/3*u[1])
end

function g(du, u, p, t)
    du[1] = 3/5
end

# Set the initial condition
u0 = [2.0]
# Set the time span
tspan = (0.0, sqrt(2))
# Set the number of realizations
nrealizations = 10000
# Set the number of time steps
nsteps = 10000
# Set the time step size
dt = (tspan[2] - tspan[1])/nsteps
# Set the number of time steps to reach the time t = 2
ntsteps = 2/dt

function WienerProcess(nsteps, nrealizations)
    # Initialize the Wiener process
    W = zeros(nsteps+1, nrealizations)
    # Initialize the random number generator
    rng = MersenneTwister(1234)
    # Generate the Wiener process
    for j in 1:nrealizations
        for i in 1:nsteps
            W[i+1, j] = W[i, j] + sqrt(dt)*randn(rng)
        end
    end
    return W
end

# Generate the Wiener process
W = WienerProcess(nsteps, nrealizations)

# Initialize the solution
X = zeros(nsteps+1, nrealizations)
# Initialize the solution for the expected value
EX = zeros(nsteps+1)
# Initialize the solution for the expected value of X^0.6
EX06 = zeros(nsteps+1)
# Initialize the solution for the probability
PX = zeros(nsteps+1)

# Generate the realizations
for j in 1:nrealizations
    # Initialize the solution
    X[1, j] = u0[1]
    # Generate the realizations
    for i in 1:nsteps
        X[i+1, j] = X[i, j] + f([0.0], [X[i, j]], [], 0.0)[1]*dt + g([0.0], [X[i, j]], [], 0.0)[1]*sqrt(dt)*randn()
    end
end

# Compute the expected value
for i in 1:nsteps+1
    EX[i] = mean(X[i, :])
    # EX06[i] = mean(X[i, :].^0.6) make sure it's a real number
    if mean(X[i, :].^0.6) < Inf
        EX06[i] = mean(X[i, :].^0.6)
    else
        EX06[i] = 0
    end
    PX[i] = sum(X[i, :] .> 2)/nrealizations
end


plot!(0:dt:tspan[2], EX06, label=L"$E[X^{0.6}]$", xlabel=L"$t$", ylabel=L"$E[X^{0.6}]$", title=L"Expected Value of $X^{0.6}$", legend=:topleft)

# Save the plot
savefig("imgs/question3_plot.png")