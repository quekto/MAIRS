# ==============
# Version 1.0.0
# Created by Ian Fitch Mochida at University of Tokyo
# Calculations based on Mastrapa et al, 2019 and Hudgins et al, 1992
# ==============
using Optim
using CSV
using DataFrames
using Plots
using Printf
using Dates
using Statistics
using DelimitedFiles
using PyCall
using Conda
# Conda.add("scipy")
# Constants
const a::Float64 = 0.29
# const d::Float64 = 0.9e-8 # m
const d_cent::Float64 = 0.9e-8 * 100 # cm
const H::Float64 = 0.334 # Si
const n_b::Float64 = 1.32

const file_path::String = "/Users/ianfitch/Library/CloudStorage/GoogleDrive-ian-fitch-mochida@g.ecc.u-tokyo.ac.jp/My Drive/2_Lab/MAIRS/2024/0403_10K_H2O_32min/0403IPOPow.csv"
#  "/Users/ianfitch/Library/CloudStorage/GoogleDrive-ian-fitch-mochida@g.ecc.u-tokyo.ac.jp/My Drive/2_Lab/MAIRS/2024/0922_20K_H2O_32min/0922IPOPow.csv"

# "/Users/ianfitch/Library/CloudStorage/GoogleDrive-ian-fitch-mochida@g.ecc.u-tokyo.ac.jp/My Drive/2_Lab/MAIRS/2024/0729_20K_H2O_32min/0729IPOPow.csv"
 
# "/Users/ianfitch/Library/CloudStorage/GoogleDrive-ian-fitch-mochida@g.ecc.u-tokyo.ac.jp/My Drive/2_Lab/MAIRS/2024/0725_20K_H2O_32min/0725IPOPow.csv"
# Define the wavenumber range for plotting
const lower_limit::Float64 = 3000
const upper_limit::Float64 = 4000

# Function to load CSV data
function load_data_mairs(filepath::String)
    data = CSV.read(filepath, DataFrame)
    nu_ip = data[!, 1]    # Wavenumber (cm^-1) for ip
    A_ip = data[!, 2]    # Absorbance for ip
    nu_op = data[!, 3]    # Wavenumber (cm^-1) for op
    A_op = data[!, 4]    # Absorbance for op
    return nu_ip, A_ip, nu_op, A_op
end

function load_data(filepath::String)
    data = CSV.read(filepath, DataFrame)
    nu::Vector{Float64} = data[!, 1]
    A::Vector{Float64}  = data[!, 2]
    return nu, A
end

# Function to read the external data file for "Amorphous 15"
function read_amorphous_data(file_path::String)
    # Initialize empty arrays for wavelength, n, and k
    wavelength_μm = Float64[]
    n_external = Float64[]
    k_external = Float64[]

    # Open and read the file
    file = open(file_path, "r")
    for line in eachline(file)
        if startswith(line, "Amorphous    15")  # Filter only "Amorphous 15" rows
            parts = split(line)
            push!(wavelength_μm, parse(Float64, parts[3]))  # Column 3: Wavelength in micrometers
            push!(n_external, parse(Float64, parts[4]))      # Column 4: n
            push!(k_external, parse(Float64, parts[5]))      # Column 5: k
        end
    end
    close(file)

    return wavelength_μm, n_external, k_external
end

# Function to calculate k based on the MAIRS method
function calculate_k_mairs(nu::Vector{Float64}, A_ip::Vector{Float64}, A_op::Vector{Float64}, n::Vector{Float64})
    k    = [log(10) / (16 * π * a * n[i] * d_cent * nu[i]) * (2 * A_ip[i] + A_op[i]) / 3 for i in 1:length(nu)]
    k_xy = [A_ip[i] * log(10) / (16 * π * d_cent * a * n[i] * nu[i]) for i in 1:length(nu)] # Hasegawa 2024, eq(11)
    k_z  = [A_op[i] * log(10) / (16 * π * d_cent * a * n[i] * nu[i]) for i in 1:length(nu)] # Hasegawa 2024, eq(12)
    return k, k_xy, k_z
end

function calculate_k_mairs_num(nu::Vector{Float64}, A_ip::Vector{Float64}, A_op::Vector{Float64}, n::Vector{Float64})
    # k    = [log(10) / (16 * π * a * n[i] * d_cent * nu[i]) * (2 * A_ip[i] + A_op[i]) / 3 for i in 1:length(nu)]
    k_xy::Vector{Float64} = [A_ip[i] * log(10) / (16 * π * d_cent * a * n[i] * nu[i]) for i in 1:length(nu)] # Hasegawa 2024, eq(11)
    k_z = Vector{Float64}(undef, length(nu))

    for i in eachindex(nu)
        function objective(k) # Value to minimize
            return abs(A_op[i] - (16 * pi * d_cent * a * nu[i]) / log(10) * n[i] * k[1] / (n[i] * n[i] + k[1] * k[1]) * (n[i] * n[i] + k[1] * k[1])) # 1.2^4 = 2.0736 
        end

        function grad(k::Vector)
            2 * (A_op[i] - (16 * pi * d_cent * a * nu[i]) / log(10) * n[i] * 1.20^4 * k[1] / (n[i]^2 + k^2)^2) * (16 * pi * d_cent * a * nu[i] * n[i] * 1.20^4 / log(10)) * (6 * k[1] - (n^2 + k[1]^2))/(n^2 + k[1]^2)^3
        end
        # optimizer_z = optimize(objective,  0, 1)
        # optimizer_z = optimize(objective_vec, grad,  [1.], method=NelderMead())
        # k_z_optimal = Optim.minimizer(optimizer_z)[1]

        scipy_optimize = pyimport("scipy.optimize").minimize
        result = scipy_optimize(objective, 0.5, method="L-BFGS-B", bounds=[(0,1)])
        k_z_optimal = result["x"][1]

        if i%100==0 println(i) end

        if i in 1140:1150
            # println(optimizer_z)
        end

        k_z[i] = k_z_optimal
    end
    k = (2 * k_xy + k_z) / 3
    return k, k_xy, k_z
end

function calculate_k_tr(nu::Vector{Float64}, A::Vector{Float64}, n::Vector{Float64})
    k = [log(10) * (n_1 + n_3) / (16 * pi * d_2 * n[i] * nu[i]) for i in 1:length(nu)]
end

# Function to calculate n using Kramers-Kronig relation
function calculate_n(nu::Vector{Float64}, k::Vector{Float64})
    n_new = Vector{Float64}(undef, length(nu))
    for i in 1:length(nu)
        nu_i = nu[i]
        integral = sum([4 * π * nu_j * k_j / (nu_j^2 - nu_i^2) for (nu_j, k_j) in zip(nu, k) if nu_j ≠ nu_i]) * (nu[2] - nu[1])
        n_new[i] = n_b - integral / 2pi^2 # it is n-ik in Bergren et al, 1978, but we use n+ik here. #TODO Further consideration on KK relation.
    end
    return n_new
end

# Function to calculate MAIRS spectra
function calculate_mairs_spectra(nu::Vector{Float64}, n::Vector{Float64}, k_xy::Vector{Float64}, k_z::Vector{Float64})
    A_ip_calc = [8 * π * d_cent * a * nu[i] / log(10) * 2 * n[i] * k_xy[i] for i in 1:length(nu)]
    A_op_calc = [8 * π * d_cent * a * nu[i] / log(10) * 2 * n[i] * k_z[i] for i in 1:length(nu)]
    return A_ip_calc, A_op_calc
end

# Function to check convergence
function check_convergence(A_ip_calc::Vector{Float64}, A_op_calc::Vector{Float64}, A_ip::Vector{Float64}, A_op::Vector{Float64})
    avg_diff_ip = mean(abs.(A_ip_calc .- A_ip)) / mean(A_ip)
    avg_diff_op = mean(abs.(A_op_calc .- A_op)) / mean(A_op)
    return abs((avg_diff_ip + avg_diff_op) / 2) < 0.001
end

# Main function
function mairs()
    # Load data
    nu_ip, A_ip, nu_op, A_op = load_data_mairs(file_path) # Replace with actual file path
    A_op = A_op
    # nu_ip, A_ip, nu_op, A_op = nu_ip[1:6200], A_ip[1:6200], nu_op[1:6200], A_op[1:6200]
    # n::Float64 = 1
    n_values::Vector{Float64} = ones(length(nu_ip))
    k::Vector{Float64} = ones(length(nu_ip))   # Start with k as ones

    # Create unique folder for plots
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    folder = "fig_$timestamp"
    mkpath(folder)

    # Load external data
    wavelength_μm, n_external, k_external = read_amorphous_data("apj315328t8_mrt.txt")

    # Convert wavelength in micrometers to wavenumber in cm⁻¹ for comparison
    nu_external = 1.0 ./ wavelength_μm * 1e4  # cm⁻¹

    # Filter external data based on the specified range
    external_filtered_indices = (nu_external .>= lower_limit) .& (nu_external .<= upper_limit)
    nu_external_filtered = nu_external[external_filtered_indices]
    n_external_filtered = n_external[external_filtered_indices]
    k_external_filtered = k_external[external_filtered_indices]

    # Iterative calculation loop
    iteration = 1
    while true
        println("Iteration $iteration")

        # Calculate k
        k, k_xy, k_z = calculate_k_mairs_num(nu_ip, A_ip, A_op, n_values)

        # Update n using KK relation
        n_values = calculate_n(nu_ip, k)

        # Calculate MAIRS spectra
        A_ip_calc, A_op_calc = calculate_mairs_spectra(nu_ip, n_values, k_xy, k_z)

        # Plot and save
        plot(nu_ip, A_ip, label="A_ip (Input)", color=:blue, title="MAIRS - $(basename(file_path))")
        plot!(nu_ip, A_ip_calc, label="A_ip (Calculated)", color=:red, linestyle=:dash)
        plot!(nu_op, A_op, label="A_op (Input)", color=:green)
        plot!(nu_op, A_op_calc, label="A_op (Calculated)", color=:purple, linestyle=:dash)
        xflip!(true)
        savefig(joinpath(folder, "MAIRS_iteration_$iteration.png"))

        # Plot n on the left y-axis
        # Filter the data based on the wavenumber range
        filtered_indices = (nu_ip .>= lower_limit) .& (nu_ip .<= upper_limit)
        nu_filtered = nu_ip[filtered_indices]
        n_filtered = n_values[filtered_indices]
        k_filtered = k[filtered_indices]

        # plt = plot()
        plot(nu_filtered, k_filtered, label="k", color=:red, ylabel="Extinction Coefficient (k)", title="(n) and (k) - $(basename(file_path)), d=$(d_cent/100) m")
        right_axis = twinx()
        plot!(right_axis, nu_filtered, n_filtered, label="n", color=:blue, xlabel="Wavenumber (cm⁻¹)", ylabel="Refractive Index (n)") # Plot k on the right y-axis
        # Plot external n and k
        plot!(nu_external_filtered, k_external_filtered, label="k (External)", color=:orange, linestyle=:dash)
        plot!(right_axis, nu_external_filtered, n_external_filtered, label="n (External)", color=:cyan, linestyle=:dash)
        # Reverse the x-axis
        xflip!(true)
        plot!(legend=:topright)

        # Save the figure to the specified folder
        savefig( joinpath(folder, "n_k_iteration_$iteration.png"))

        # Check for convergence
        if check_convergence(A_ip_calc, A_op_calc, A_ip, A_op)
            println("Convergence achieved. Process completed.")
            break
        end

        iteration += 1
    end
end

function Tr()
    nu, A = load_data(file_path)

    n_values::Vector{Float64} = ones(length(nu_ip))
    k::Vector{Float64} = ones(length(nu_ip))   # Start with k as ones

    # Create unique folder for plots
    timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    folder = "fig_$timestamp"
    mkpath(folder)

    iteration = 1
    while true
        println("Iteration $iteration")
        k = calculate_k_tr(nu)

        n_values = calculate_k
    end

end

# Run the main function
mairs()


# Numerical Calculation
function reverse(k)
    println(A_op[5205] * log(10) / (16pi * d_cent * a * nu_op[5205]) - n * k / (n^2 + k^2)^2)
    println(16 * pi * d_cent * a / log(10) * (n * k / (n^2 + k^2)^2))
end



function temp()
    nu_ip, A_ip, nu_op, A_op = load_data_mairs(file_path) # Replace with actual file path
    n = 1.1
    function objective(k)
        (A_op[5205] * log(10) / (16pi * d_cent * a * nu_op[5205]) - n * k / (n^2 + k^2)^2 ) ^2
    end
    function objective_2(k)
        (A_op[5205] * log(10) / (16 * π * d_cent * a * n * nu_op[5205]) - k)^2
    end

    result = optimize(objective, 0.25, 1)

    k_optimal = Optim.minimizer(result)
    
    # println(cal_obj())
    println("Optimal k value: ", k_optimal)
    reverse(k_optimal)
    function f(k_optimal)
        return A_op[5205] * log(10) / (16pi * d_cent * a * nu_op[5205]) - n * k_optimal / (n^2 + k_optimal^2)^2
    end

end

# temp()