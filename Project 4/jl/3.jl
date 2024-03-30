
using Random
using Statistics
using Plots
using LaTeXStrings
using DifferentialEquations

# Consider the following SDE
# 𝑑𝑋(𝑡) = 𝜇𝑋(𝑡)𝑑𝑡 + 𝜎𝑋(𝑡)𝑑𝑊(𝑡) , 𝑋(0) = 3, 𝜇 = 2 , 𝜎 = 0.10
# Where 𝑡 ∈ [0,1]

# Parameters
μ = 2
σ = 0.10
X₀ = 3
tspan = (0.0, 1.0)

# Define the SDE problem
function my_SDE(du, u, p, t)
    du[1] = μ*u[1]
    du[2] = σ*u[2]
end

# Define the noise process
function noise(du, u, p, t)
    du[1] = 0.0
    du[2] = 1.0
end

# Define the initial condition
u₀ = [X₀, 0.0]
# Define the number of points
N = 1000
# Define the number of simulations
M = 1000
# Define the time points
t = range(tspan[1], tspan[2], length=N)

function Wiener_process(N)
    Δt = t[2] - t[1]
    dW = sqrt(Δt)*randn(N)
    W = cumsum(dW)
    return W
end

# Define the Wiener process
W = Wiener_process(N)

# Define the solution
sol = solve(SDEProblem(my_SDE, noise, u₀, tspan), EM(), dt=1e-3, saveat=t)

# Plot the solution
plot(sol, idxs = 1, label = L"X(t)", xlabel = L"t", ylabel = L"X(t)", title = "SDE Solution", color = :blue, lw = 2)

# Plot the Wiener process
plot!(t, W, label = L"W(t)", color = :red, lw = 2)

# Save the plot
savefig("./imgs/3.png")

# a) Show that the Euler Maruyama method has weak order of convergence equal to one. That is
# |𝐸[𝑋1] − 𝐸[𝑋(1)]| = 𝐶Δ𝑡. Here 𝑋(1) is the exact solution at time 1 and 𝑋1 is the computed solution at time
# 1

# Define the error
error = zeros(M)
delta_t = 1e-3

function simulate()
    # Define the Wiener process
    W = Wiener_process(N)
    # Define the solution
    sol = solve(SDEProblem(my_SDE, noise, u₀, tspan), EM(), dt=delta_t, saveat=t)
    # Compute the error
    return abs(sol(1.0)[1] - 3)
end
