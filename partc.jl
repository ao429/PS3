###### Worked with Azmain ######

include("Flux.jl")
using LinearAlgebra
using CSV
S = Matrix{Float64}(CSV.read("stoichiometric_matrix.csv", header=0))

(species,fluxes) = size(S);
default_bounds_array = zeros(fluxes,2);

objective_coefficient_array = zeros(fluxes);

species_bounds_array = zeros(species,2);

k_cat = [88.1;           ## 2.1.3.3
        34.5;           ## 4.3.2.1
         249;           ## 3.5.3.1
        203;           ## 6.3.4.5
        13.7;           ## 1.14.13.39
        13.7];            ## 1.14.13.39              

E = 0.01                ## umol/gDW

c = [0.923*0.9897*1;   ## 6.3.4.5
      1;                ## 4.3.2.1
      0.1418;           ## 3.5.3.1
      0.7372*1;         ## 2.1.3.3
      1*0.9865;         ## 1.14.13.39
      1;                ## 1.14.13.39
    ];                  ## Park 

v_max = 3600/10^3*E.*k_cat.*c;  

default_bounds_array[1:length(v_max),2] .= v_max;
default_bounds_array[(length(v_max)+1):end,2] .= 10;
default_bounds_array[(15:20),1] .= -10;

objective_coefficient_array[10] = -1;

optimize = (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag)=calculate_optimal_flux_distribution(stoichiometric_matrix, default_bounds_array, species_bounds_array, objective_coefficient_array);
optimized_flux=optimize[2];
max_urea_flux=optimized_flux[10];

println("maximum urea flux = ", max_urea_flux, " mmol/gDW-hr")