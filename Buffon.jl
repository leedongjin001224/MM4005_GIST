#Dongjin Lee, GIST
#2024 Fall, course in Scientific Programming
#Last updated: 2024-11-15
# Program for solving Buffon's needle problem
# Resolved the problem occured when executed on ssh server by changing the libraries.
# Replaced the function to use quasi Monte Carlo method.

# Command "julia ver6.jl" to run this code.


# Defining parameters and variables.
log10_N = 8 # Number of array elements at one iteration.
n = 10 # Iterating number.
print("Enter the length between two neighbor lines:")
d = readline() # Length between two lines.
d = parse(Int,d)
print("Enter the length of the needle:")
l = readline() # Length of the needle.
l = parse(Int,l)


print("Loading libraries.\n")
using Statistics
using GLM
import DataFrames
using Plots
using LaTeXStrings
using QuasiMonteCarlo
using Distributions
print("Loaded libraries.\n")



# N is the size of the each one sample.
# n is the number of samples.
function solve_buffon(N, n, d, l, pi_)
    sample_means = zeros(n)
    for i in 1:n
        samples = QuasiMonteCarlo.sample(N,2,Uniform())
        x = d * (1 .- samples[1,:])
        # Twice of the distance between the center of the needle and closest line from it.
        θ = samples[2,:] * π / 2
        # Angle the needle rotated. 0 corresponds to horizontal line.

        trials = x .< l * cos.(θ)
        # 1 if needle touched the line, else 0.
        sample_means[i] = sum(trials) / N
        # Leave only sample mean.
    end

    return sum(sample_means) / n
end

function π_error_plotting(method, results, n)
    # Change plotting index and data into log scale.
    abs_errors = abs.(results .- π)[:]

    N_indices = zeros(0)
    for i in 1:log10_N
        N_indices = [
            N_indices
            ones(n) * (10 ^ i)
        ]
    end

    log10_N_indices = log10.(N_indices)
    log10_abs_errors = log10.(abs_errors)

    # Regression line below.
    data = DataFrames.DataFrame(x=log10_N_indices, y=log10_abs_errors)
    model = lm(@formula(y ~ x), data)
    coefs = coef(model)
    r_2 = r2(model)
    ste = stderror(model)



    # Plotting and saving below.

    scatter(log10_N_indices, log10_abs_errors, alpha=0.2, markersize=10, label=L"\mathrm{Trials}", title = "Monte Carlo Result for $(method)")
    plot!([0,log10_N+1], [coefs[1], coefs[1] + coefs[2] * (log10_N+1)], color="red", label="slope:$(coefs[2])+-$(2*ste[2]) R^2=$(r_2)")
    xlabel!(L"x=\log N,\mathrm{\,size\,of\,the\,sample\,in\,log-scale}",)
    ylabel!(L"y=\log |\hat{\pi}-\pi|,\mathrm{\,absolute\,error\,in\,log-scale}")
    #axislegend(ax, merge=true)
    savefig("$(method).pdf")
    print("Plot saved as a pdf file.\n")
    ### Save the plot as a pdf file because it is not displayed in cmd.
end

## Pi approximation with the Hit-and-Miss algorithm
print("Pi approximating with the Hit-and-Miss algorithm.\n")
results = zeros((n,log10_N))

for l in 1:log10_N
    global π_sample_means = zeros(n)
    N = 10^l
    for i in 1:n
        samples = QuasiMonteCarlo.sample(N, 2, Uniform())
        r = samples[1,:] .^ 2 + samples[2,:] .^ 2
        π_sample_means[i] = sum(r .< 1) * 4 / N
    end
    results[:,l] = π_sample_means
end

π_h = sum(π_sample_means) / n

π_error_plotting("Hit-and-Miss_approximation", results, n)

print("Hit-and-Miss approximation to pi: ", π_h, "\n")
print("Solution of the Buffon's needle problem with approximated pi: ", solve_buffon(10^log10_N, n, d, l, π_h), "\n")



## Pi approximation with the Wallis product
print("Pi approximating with the Wallis product.\n")
π_w = 2
pivot = 1
results = zeros((1,log10_N))
for i in 1:10^log10_N
    global π_w *= 4 * i^2
    global π_w /= 4 * i^2 - 1
    if 10^pivot == i
        results[1,pivot] = π_w
        global pivot +=1
    end
end

π_error_plotting("Wallis_product", results, 1)

print("Wallis product approximation to pi: ", π_w, "\n")
print("Solution of the Buffon's needle problem with approximated pi: ", solve_buffon(10^log10_N, n, d, l, π_w), "\n")


## Bibliographies   
#[1] HWANG, Chi-Ok, et al. Buffon’s Needle Algorithm to Estimate π. *Applied Mathematics*, 2017, 8.3: 275-279.   
#[2] WALLIS, John; STEDALL, Jacqueline A. *The arithmetic of infinitesimals*. New York: Springer, 2004.
