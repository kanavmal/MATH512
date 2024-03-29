# Estimate the following expected value: E(W_3^2 + sin(W_3) + 2e^{W_3})
# Where W_t is a standard Wiener process, that is the drift parameter is zero and the Variance parameter sigma^2 = 1.

# The expected value of a function of a Wiener process is given by the following formula:
# E(f(W_t)) = \int_{-\infty}^{\infty} f(x) * \frac{1}{\sqrt{2\pi t}} * e^{-\frac{x^2}{2t}} dx
# Where t is the time at which we are evaluating the function f, and x is the value of the Wiener process at time t.

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

function expected_value(n)
    # Generate n samples from a standard Wiener process
    w = wiener_process(n)
    # Calculate the expected value of the function f(W_3) = W_3^2 + sin(W_3) + 2e^{W_3}
    f = w.^2 .+ sin.(w) .+ 2 .* exp.(w)
    return mean(f)
end

n = 200000
expected_value(n)

println("Estimated: ", expected_value(n))

# plot the histogram of the function f(W_3) = W_3^2 + sin(W_3) + 2e^{W_3}
w = wiener_process(n)
f = w.^2 .+ sin.(w) .+ 2 .* exp.(w)
# scale the x axis from 0, 50
# histogram(f, bins=50, label=L"$f(W_3)$", xlabel=L"$f(W_3)$", ylabel="Frequency", title=L"Histogram of $f(W_3)$")
histogram(f, bins=50, label=L"$f(W_3)$", xlabel=L"$f(W_3)$", ylabel="Frequency", title=L"Histogram of $f(W_3)$", xlims=(0, 300))
savefig("./imgs/question1_histogram.png")

# We plot the expected value of the function as it varies with the
#     number of samples in Figure 2.
#     Figure: Convergence of the Given Expectation Over 200,000 Iterations
    
function convergence(n)
    # Calculate the expected value of the function f(W_3) = W_3^2 + sin(W_3) + 2e^{W_3}
    f = w.^2 .+ sin.(w) .+ 2 .* exp.(w)
    # Calculate the expected value of the function f(W_3) = W_3^2 + sin(W_3) + 2e^{W_3}
    expected_values = [mean(f[1:i]) for i in 1:n]
    return expected_values
end

n = 200000
expected_values = convergence(n)
plot(expected_values, label="Expected Value", xlabel="Number of Samples", ylabel="Expected Value")
# add a horizontal line for the converged value
hline!([expected_value(n)], label="Converged Value", color="red")
savefig("./imgs/question1_convergence.png")


# find the contributing term that is explaining the most variance
# W_t^2, sin(W_t), 2e^{W_t}, calculate each function's mean and variance
w = wiener_process(n)
w_squared = w.^2
w_sin = sin.(w)
w_exp = 2 .* exp.(w)

println("Mean of W_t^2: ", mean(w_squared))
println("Variance of W_t^2: ", var(w_squared))

println("Mean of sin(W_t): ", mean(w_sin))
println("Variance of sin(W_t): ", var(w_sin))

println("Mean of 2e^{W_t}: ", mean(w_exp))
println("Variance of 2e^{W_t}: ", var(w_exp))