# Let S𝑡 be a Geometric Brownian Motion process: S_t = S_0e^{sigma W_t + (r-sigma^2/2)t} 
# where r = 0.05, sigma = 0.2, S_0 = 90, and W_t is a standard Wiener process. Estimate
# E[S_3].

using Random
using Statistics
using Plots
using LaTeXStrings

function wiener_process(n)
    # Generate n samples from a standard normal distribution
    z = randn(n)
    # Calculate the Wiener process using the formula W_t = sqrt(t) * Z
    w = sqrt.(3) .* z
    return w
end

function geometric_brownian_motion(n)
    # Generate n samples from a standard Wiener process
    w = wiener_process(n)
    # Calculate the Geometric Brownian Motion process using the formula S_t = S_0e^{sigma W_t + (r-sigma^2/2)t}
    r = 0.05
    sigma = 0.2
    S_0 = 90
    t = 3
    S = S_0 .* exp.(sigma .* w .+ (r - sigma^2 / 2) .* t)
    return S
end

function expected_value(n)
    # Generate n samples from a Geometric Brownian Motion process
    S = geometric_brownian_motion(n)
    # Calculate the expected value of the Geometric Brownian Motion process
    return mean(S)
end

n = 200000
expected_value(n)

println("Estimated: ", expected_value(n))
# calculate the real value of E[S_3]
r = 0.05
sigma = 0.2
S_0 = 90
t = 3
real_value = S_0 * exp(r * t)
println("Real Value: ", real_value)