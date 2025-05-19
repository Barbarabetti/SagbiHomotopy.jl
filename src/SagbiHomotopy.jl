module SagbiHomotopy

using Oscar
using HomotopyContinuation

export wDeg, 
get_coeffs_exponents, 
initial_form, 
SagbiCriterion, 
#extractWeightVectors,
weightVectorsRealizingSAGBI,
HC_to_oscar,
get_weight,
sagbi_homotopy,
possible_sagbi,
get_leading_monomial,
weight_deformation_for_poly,
degree_map,
degree_monomial_map,
isolated_nsols_baselocus,
get_base_locus,
get_sagbi_grassmannian


include("Utils.jl")
include("SagbiDetection.jl")
include("Homotopy.jl")

end # module SagbiHomotopy