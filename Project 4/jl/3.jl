
using Random
using Statistics
using Plots
using LaTeXStrings
using DifferentialEquations
using Polynomials

# Consider the following SDE
# 𝑑𝑋(𝑡) = 𝜇𝑋(𝑡)𝑑𝑡 + 𝜎𝑋(𝑡)𝑑𝑊(𝑡) , 𝑋(0) = 3, 𝜇 = 2 , 𝜎 = 0.10
# Where 𝑡 ∈ [0,1]

# a) Show that the Euler Maruyama method has weak order of convergence equal to one. That is
# |𝐸[𝑋1] − 𝐸[𝑋(1)]| = 𝐶Δ𝑡. Here 𝑋(1) is the exact solution at time 1 and 𝑋1 is the computed solution at time
# 1

x0 = 3
μ = 2
σ = 0.10
t = 1
dt = 0.01
n = Int(t / dt)
sample_dts = 2 .^ (5:-1:1) .* dt
number_of_samples = 10000

function euler_maruyama(a, dt, dw)
    x = a
    return x + μ * x * dt + σ * x * dw
end

function exact_solution(x0, μ, σ, t, dw)
    return x0 * exp((μ - 0.5 * σ^2) * t + σ * dw)
end

function compute_error(x0, μ, σ, t, dt, n, number_of_samples)
    errors = zeros(number_of_samples)
    for i in 1:number_of_samples
        dw = sqrt(dt) * randn()
        x = x0
        for j in 1:n
            x = euler_maruyama(x, dt, dw)
            dw = sqrt(dt) * randn()
        end
        errors[i] = abs(x - exact_solution(x0, μ, σ, t, dw))
    end
    return mean(errors)
end

errors = [compute_error(x0, μ, σ, t, dt, n, number_of_samples) for dt in sample_dts]

