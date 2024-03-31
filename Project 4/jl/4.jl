# Consider the following SDE:
# 𝑑𝑋(𝑡) = 𝜇𝑋(𝑡)𝑑𝑡 + 𝜎𝑋(𝑡)𝑑𝑊(𝑡) , 𝑋(0) = 3, 𝜇 = 2 , 𝜎 = 0.10
# a) Simulate (over the interval [0,20]) this stochastic process using an implicit method of the form
# 𝑋𝑛+1 = 𝑋𝑛 + (1 − 𝜃)𝛥𝑡𝑓(𝑋𝑛) + 𝜃𝛥𝑡𝑓(𝑋𝑛+1) + √𝛥𝑡𝛼𝑛𝑔(𝑋𝑛)

using Plots
using Random
using Statistics
using LaTeXStrings
using DifferentialEquations

function f(x, mu)
    return mu * x
end

function g(x, sigma)
    return sigma * x
end

function implicit_method(mu, sigma, x0, t, dt, theta)
    n = Int(t / dt)
    x = zeros(n + 1)
    x[1] = x0
    for i in 1:n
        dw = sqrt(dt) * randn()
        x[i + 1] = x[i] + (1 - theta) * dt * f(x[i], mu) + theta * dt * f(x[i + 1], mu) + sqrt(dt) * g(x[i], sigma) * dw
    end
    return x
end


mu = 2
sigma = 0.10
x0 = 3
t = 20
dt = 0.01
theta = 0.5

x = implicit_method(mu, sigma, x0, t, dt, theta)
plot(0:dt:t, x, label=L"X(t)", xlabel="Time", ylabel="Value", title=L"Stochastic process $X(t)$", legend=:topleft)
savefig("./imgs/4a.png")


# b) Compare with the analytical solution.

function analytical_solution(mu, sigma, x0, t, dt)
    n = Int(t / dt)
    x = zeros(n + 1)
    x[1] = x0
    for i in 1:n
        dw = sqrt(dt) * randn()
        x[i + 1] = x0 * exp((mu - 0.5 * sigma^2) * i * dt + sigma * dw)
    end
    return x
end

x_analytical = analytical_solution(mu, sigma, x0, t, dt)
plot!(0:dt:t, x_analytical, label=L"$X(t)$ (Analytical)", xlabel="Time", ylabel="Value", title=L"Stochastic process $X(t)$", legend=:topleft)
savefig("./imgs/4b.png")

# c) For what values of 𝜇 𝑎𝑛𝑑 𝜎 is the SDE mean-square stable.

# The SDE is mean-square stable if the following condition is satisfied:
# It is called mean-square stable if, for every ε > 0, there exists δ > 0 such that
# E[|X(t)|^2] ≤ ε for all t ≥ 0 whenever E[|X(0)|^2] ≤ δ.
# Let's test this condition


# Test for mean-square stability
function test_mean_square_stability(mu, sigma, x0, t, dt, theta)
    n = Int(t / dt)
    x = zeros(n + 1)
    x[1] = x0
    for i in 1:n
        dw = sqrt(dt) * randn()
        x[i + 1] = x[i] + (1 - theta) * dt * f(x[i], mu) + theta * dt * f(x[i + 1], mu) + sqrt(dt) * g(x[i], sigma) * dw
    end
    return mean(x.^2)
end

# Test for different values of mu and sigma
mu_values = 0:0.1:5
sigma_values = 0.01:0.01:0.5
stable_mu_sigma = []
for mu in mu_values
    for sigma in sigma_values
        x = test_mean_square_stability(mu, sigma, x0, t, dt, theta)
        if x < 1e6
            push!(stable_mu_sigma, (mu, sigma, x))
        end
    end
end

# print the range for which mu and sigma are mean-square stable
range_mu = (minimum([x[1] for x in stable_mu_sigma]), maximum([x[1] for x in stable_mu_sigma]))
range_sigma = (minimum([x[2] for x in stable_mu_sigma]), maximum([x[2] for x in stable_mu_sigma]))
println("Mean-square stable range for mu: ", range_mu)
println("Mean-square stable range for sigma: ", range_sigma)

# d) For what values of 𝜃 is the implicit method mean-square stable.

# The implicit method is mean-square stable if the following condition is satisfied:
# It is called mean-square stable if, for every ε > 0, there exists δ > 0 such that
# E[|X(t)|^2] ≤ ε for all t ≥ 0 whenever E[|X(0)|^2] ≤ δ.
# Let's test this condition

# Test for mean-square stability
function test_mean_square_stability(mu, sigma, x0, t, dt, theta)
    n = Int(t / dt)
    x = zeros(n + 1)
    x[1] = x0
    for i in 1:n
        dw = sqrt(dt) * randn()
        x[i + 1] = x[i] + (1 - theta) * dt * f(x[i], mu) + theta * dt * f(x[i + 1], mu) + sqrt(dt) * g(x[i], sigma) * dw
    end
    return mean(x.^2)
end

# Test for different values of theta
theta_values = 0:0.1:1
stable_theta = []
for theta in theta_values
    x = test_mean_square_stability(mu, sigma, x0, t, dt, theta)
    if x < 1e6
        push!(stable_theta, (theta, x))
    end
end

# print the range for which theta is mean-square stable
range_theta = (minimum([x[1] for x in stable_theta]), maximum([x[1] for x in stable_theta]))
println("Mean-square stable range for theta: ", range_theta)

# e) For what values of 𝜇 𝑎𝑛𝑑 𝜎 is the SDE asymptotically stable.

# The SDE is asymptotically stable if the following condition is satisfied:
# It is called asymptotically stable if, for every ε > 0, there exists δ > 0 such that
# lim t→∞ E[|X(t)|^2] ≤ ε whenever E[|X(0)|^2] ≤ δ.
# Let's test this condition

# Test for asymptotic stability
function test_asymptotic_stability(mu, sigma, x0, t, dt, theta)
    n = Int(t / dt)
    x = zeros(n + 1)
    x[1] = x0
    for i in 1:n
        dw = sqrt(dt) * randn()
        x[i + 1] = x[i] + (1 - theta) * dt * f(x[i], mu) + theta * dt * f(x[i + 1], mu) + sqrt(dt) * g(x[i], sigma) * dw
    end
    return mean(x.^2)
end

# Test for different values of mu and sigma
mu_values = 0:0.1:5
sigma_values = 0.01:0.01:0.5
stable_mu_sigma = []
for mu in mu_values
    for sigma in sigma_values
        x = test_asymptotic_stability(mu, sigma, x0, t, dt, theta)
        if x < 1e6
            push!(stable_mu_sigma, (mu, sigma, x))
        end
    end
end

# print the range for which mu and sigma are asymptotically stable
range_mu = (minimum([x[1] for x in stable_mu_sigma]), maximum([x[1] for x in stable_mu_sigma]))
range_sigma = (minimum([x[2] for x in stable_mu_sigma]), maximum([x[2] for x in stable_mu_sigma]))
println("Asymptotically stable range for mu: ", range_mu)
println("Asymptotically stable range for sigma: ", range_sigma)

# f) For what values of 𝜃 is the Implicit method asymptotically stable
 
# The implicit method is asymptotically stable if the following condition is satisfied:
# It is called asymptotically stable if, for every ε > 0, there exists δ > 0 such that
# lim t→∞ E[|X(t)|^2] ≤ ε whenever E[|X(0)|^2] ≤ δ.
# Let's test this condition

# Test for asymptotic stability
function test_asymptotic_stability(mu, sigma, x0, t, dt, theta)
    n = Int(t / dt)
    x = zeros(n + 1)
    x[1] = x0
    for i in 1:n
        dw = sqrt(dt) * randn()
        x[i + 1] = x[i] + (1 - theta) * dt * f(x[i], mu) + theta * dt * f(x[i + 1], mu) + sqrt(dt) * g(x[i], sigma) * dw
    end
    return mean(x.^2)
end

# Test for different values of theta
theta_values = 0:0.1:1
stable_theta = []
for theta in theta_values
    x = test_asymptotic_stability(mu, sigma, x0, t, dt, theta)
    if x < 1e6
        push!(stable_theta, (theta, x))
    end
end

# print the range for which theta is asymptotically stable
range_theta = (minimum([x[1] for x in stable_theta]), maximum([x[1] for x in stable_theta]))
println("Asymptotically stable range for theta: ", range_theta)