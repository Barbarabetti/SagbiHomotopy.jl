module SagbiHomotopy

import Oscar
import HomotopyContinuation
using HomotopyContinuation: solutions
using HomotopyContinuation.ModelKit: @unique_var, @var, Expression, System, Variable, exponents_coefficients, support_coefficients
using LinearAlgebra: det, dot, norm, I
using Combinatorics: combinations
using MultivariatePolynomials: subs, variables
using Oscar.Orderings: wdeglex
using AbstractAlgebra: gens, hom, ideal, kernel, matrix, polynomial_ring, quo, vars
using AbstractAlgebra.Generic: MPolyBuildCtx, finish, poly, push_term!
using Nemo: QQ, QQMPolyRing, QQMPolyRingElem
using Oscar: absolute_primary_decomposition, exponents, maximal_cones, normal_fan
using Hecke: support

export wDeg,
    get_coeffs_exponents,
    initial_form,
    SagbiCriterion,
    #extractWeightVectors,
    weightVectorsRealizingSAGBI,
    HC_to_oscar,
    get_weight,
    get_sagbi_grassmannian,
    sagbi_homotopy,
    possible_sagbi,
    get_leading_monomial,
    weight_deformation_for_poly,
    degree_map,
    degree_monomial_map,
    isolated_nsols_baselocus,
    get_base_locus


include("Utils.jl")
include("SagbiDetection.jl")
include("Homotopy.jl")

end # module SagbiHomotopy
