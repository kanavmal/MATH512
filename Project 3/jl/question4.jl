# Consider the following SDE:
# dX_t = aX_tdt + bX_tdW_t, X_0 = 100, a = 0.07, b = 0.12

# a) Simulate this stochastic process using the discretization schemes of Euler-Maruyama
# b) Compare with the analytical solution.

using Random
using Statistics
using Plots
using LaTeXStrings
using DifferentialEquations

# Define the SDE
function f(u, p, t)
    return 0.07*u
end

function g(u, p, t)
    return 0.12*u
end

# Euler-Maruyama discretization scheme
function EulerMaruyama(u0, tspan, nsteps, nrealizations)
    # Set the time step size
    dt = (tspan[2] - tspan[1])/nsteps
    # Initialize the solution
    u = zeros(nsteps+1, nrealizations)
    # Initialize the random number generator
    rng = MersenneTwister(1234)
    # Generate the realizations
    for j in 1:nrealizations
        # Set the initial condition
        u[1, j] = u0
        for i in 1:nsteps
            u[i+1, j] = u[i, j] + f(u[i, j], [], i*dt)*dt + g(u[i, j], [], i*dt)*sqrt(dt)*randn(rng)
        end
    end
    return u
end


# Set the initial condition
u0 = 100.0
# Set the time span
tspan = (0.0, 10.0)
# Set the number of time steps
nsteps = 1000
# Set the number of realizations
nrealizations = 1000

# Generate the realizations
u = EulerMaruyama(u0, tspan, nsteps, nrealizations)

# Plot the realizations
plot(0:nsteps, u, label="", xlabel="Time", ylabel="X", title="Realizations of the Stochastic Process", color=:black)
savefig("./imgs/question4.png")
