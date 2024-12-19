##### Dongjin Lee, GIST, course on Scientific Programming
##### This is a program for plotting result of Brownian motion simulation to find solve the boundary problem of the Laplace equation.



println("Plotting code is now working...")

using StatsBase, Plots, LaTeXStrings
using Statistics
using GLM
import DataFrames
using DelimitedFiles

### Read parameters
info = readdlm("./results/info.txt")
d = info[1]
ϵ = info[2]
initial = convert(Int, info[3])
log10_N = convert(Int, info[4])
N = 10^log10_N
bin_num = convert(Int, info[5])
n = convert(Int, info[6])
bin_d_ratio = info[7]



### Function for calculating radial distance from horizontal position.
function distance_radial(horizontal, d=d, input_radial=false)
    return (horizontal^2 + d^2)^0.5
end

### Function for calculating origin cummulative density function
function cdf_inverse(edges, d=d, symmetry=false)
    #arctan = atan.(edges/d)
    if symmetry
        return d * (1 .- 1 ./ (d^2 .+ edges.^2).^0.5)
    else
        return 0 #arctan / π + 0.5
    end
end

cut = bin_d_ratio * d * (bin_num-1)
ordinary_bin_width = bin_d_ratio * d

# Point where abs_errors would be calculated.
error_checking_bin_num = bin_num ÷ 3
x_ = (error_checking_bin_num-1) * cut / (bin_num-1) + 0.5 * ordinary_bin_width
r_ = distance_radial(x_)
pdf_ = d ./ r_.^2


# Walking many times.
### Now replaced by loading results from c++.

#println("Random walking...")
abs_errors = zeros((n, log10_N-initial+1))
for i in initial:log10_N
    for j in 1:n
        local filename = "./results/$(j-1)_$(i).txt"
        local bin_result = readdlm(filename)
        abs_errors[j,i-initial+1] = bin_result[error_checking_bin_num]/(10^i)/ordinary_bin_width.*r_./x_
    end
end
result = readdlm("./results/0_$(log10_N).txt")[:,1]

# Errors
abs_errors = abs.(abs_errors .- pdf_)[:]
log10_abs_errors = log10.(abs_errors)
log10_N_indices = zeros(0)
for i in initial:log10_N
    global log10_N_indices = [
        log10_N_indices
        ones(n) * i
    ]
end

println("Now plots would be saved.\nThe code might raise error- latex: command not found\nThough, plots would be saved safely.")
# Regression line.
data = DataFrames.DataFrame(x=log10_N_indices, y=log10_abs_errors)
model = lm(@formula(y ~ x), data)
coefs = coef(model)
r_2 = r2(model)
ste = stderror(model)

# Plotting and saving below.

scatter(log10_N_indices, log10_abs_errors, alpha=0.2, markersize=10, label=L"\mathrm{Trials}", title = "Monte Carlo Result for PDF at  Radial Distance $(r_)")
plot!([initial, log10_N+1], [coefs[1] + coefs[2] * initial, coefs[1] + coefs[2] * (log10_N+1)], color="red", label="slope:$(coefs[2])+-$(2*ste[2]) R^2=$(r_2)")
xlabel!(L"x=\log N,\mathrm{\,size\,of\,the\,sample\,in\,log-scale}",)
ylabel!(L"y=\log |\hat{p}-p|,\mathrm{\,absolute\,error\,in\,log-scale}")
savefig("$(bin_num)_Errors.pdf")
println("Saved plot: $(bin_num)_Errors.pdf.")



### Binning results in a dictionary
bin = Dict()
"""
bin["angular"] = Dict(
    true=> binning_angular(result, d, bin_num, true),
    false=> binning_angular(result, d, bin_num, false)
)
"""
bin["ordinary"] = Dict(
    true=> result#binning(result, d, bin_num, true),
    #false=> binning(result, d, bin_num, false)
)
edges = vcat(0:bin_num-1, Inf) * cut / (bin_num-1) .+ 0.5 * ordinary_bin_width

println("Print the index of empty bin if exists.")
#println(findall(x->x==0, bin["angular"][false].weights))
#println(findall(x->x==0, bin["angular"][true].weights))
#println(findall(x->x==0, bin["ordinary"][false].weights))
println(findall(x->x==0, bin["ordinary"][true]))


### Plots below

#p=plot(bin["ordinary"][false].weights, label="Non symmetric", dpi=300, size=(800,600), fmt=:svg)
p=plot(bin["ordinary"][true], label="Symmetric")
xlabel!("Bin index")
ylabel!("Number of particles")
title!("Number of reached particles per each bins by ordinary binning")
#display(p)
savefig(p, "$(bin_num)_Ordinary_binning.pdf")
println("Saved plot: $(bin_num)_Ordinary_binning.pdf.")

#p=plot(bin["ordinary"][false].edges[1], bin["ordinary"][false].weights/N/ordinary_bin_width, label="Non symmetric", dpi=300, size=(800,600), fmt=:svg)
p=plot(edges, bin["ordinary"][true]/N/ordinary_bin_width, label="Symmetric")
x = edges
y = x.*d./(d^2 .+ x.^2).^1.5
plot!(x, y, label=L"xd/(d^2+x^2)^{3/2}")
xlabel!("Horizontal position (absolute for symmetric plot)")
ylabel!("Probability density")
title!("PDF plot of horizontal position by ordinary binning")
lens!([0.4,0.8],[0,0.4], inset=(1,bbox(0.2,0.18,0.7,0.7)), label="")
#display(p)
savefig(p, "$(bin_num)_Horizontal_PDF.pdf")
println("Saved plot: $(bin_num)_Horizontal_PDF.pdf.")


#p=scatter(bin["ordinary"][false].edges[1], vcat(0, cumsum(bin["ordinary"][false].weights)/N), label="Non symmetric", dpi=300, size=(800,600), alpha=0.3, markerstrokewidth=0, fmt=:svg)
p=plot(edges, vcat(0, cumsum(bin["ordinary"][true])/N), label="Symmetric", alpha=0.3, markerstrokewidth=0)
#plot!(bin["ordinary"][false].edges[1], cdf_inverse.(bin["ordinary"][false].edges[1], d), label=L"\frac{1}{\pi}\arctan(x/d)+\frac{1}{2}")
plot!(edges, cdf_inverse.(edges, d, true), label=L"1-d/\sqrt{d^2+x^2}")
xlabel!("Horizontal position (absolute for symmetric plot)")
ylabel!("Cumulative probablity")
title!("CDF plot of horizontal position by ordinary binning")
lens!([0,4],[0,0.75], inset=(1,bbox(0.3,0.18,0.6,0.6)), label="")
#display(p)
savefig(p, "$(bin_num)_Horizontal_CDF_ordinary.pdf")
println("Saved plot: $(bin_num)_Horizontal_CDF_ordinary.pdf.")

#p=scatter(distance_radial.(bin["angular"][true].edges[1]), vcat(0, cumsum(bin["angular"][true].weights)/N), label="Angular binning", dpi=300, size=(800,600), alpha=0.3, markerstrokewidth=0, fmt=:svg)
p=plot(distance_radial.(edges), vcat(0, cumsum(bin["ordinary"][true])/N), label="Ordinary binning", alpha=0.3, markerstrokewidth=0)
#x = distance_radial.(edges)
#y =  d * (1 .- 1 ./ (d^2 .+ edges.^2).^0.5)
#plot!(x, y, label=L"\arccos(d/r)/\pi")
xlabel!("Radial distance")
ylabel!("Cumulative probablity")
title!("CDF plot of radial distance")
#display(p)
savefig(p, "$(bin_num)_Radial_CDF.pdf")
println("Saved plot: $(bin_num)_Radial_CDF.pdf.")


x = distance_radial.(edges[1:bin_num])
y = d ./ x.^2
p=plot(distance_radial.(edges), bin["ordinary"][true]/N/ordinary_bin_width.*x./edges[1:bin_num], label="Ordinary binning", dpi=300, size=(800,600), alpha=0.5, fmt=:svg)
plot!(x, y, label=L"d/r^2")
xlabel!("Radial distance")
ylabel!("Probability density")
title!("PDF plot of radial distance by ordinary binning")
vline!([d], linestyle=:dash, label=L"d")
lens!([0.75,1.25],[0.6,1.1], inset=(1,bbox(0.2,0.17,0.7,0.6)), label="")
#display(p)
savefig(p, "$(bin_num)_Radial_PDF.pdf")
println("Saved plot: $(bin_num)_Radial_PDF.pdf.")
