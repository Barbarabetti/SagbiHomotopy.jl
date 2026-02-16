function wDeg(monom_exponent::Vector{Int}, w::Vector{Int})
    return dot(monom_exponent, w)
end

# --- helpers for nested parametrizations ---
function normalize_param_groups(parametrization::AbstractVector)
    if isempty(parametrization)
        return Vector{Vector{Expression}}()
    elseif parametrization isa Vector{Expression}
        return [parametrization]
    elseif parametrization isa Vector{Vector{Expression}}
        return parametrization
    elseif parametrization isa Vector{Vector{Vector{Expression}}}
        return vcat(parametrization...)
    else
        if all(x -> x isa Expression, parametrization)
            return [Expression[x for x in parametrization]]
        elseif all(x -> x isa AbstractVector, parametrization)
            groups = Vector{Vector{Expression}}()
            for item in parametrization
                append!(groups, normalize_param_groups(item))
            end
            return groups
        else
            return Vector{Vector{Expression}}()
        end
    end
end

function flatten_expressions(x)
    if x isa Expression
        return Expression[x]
    elseif x isa AbstractVector
        return vcat(map(flatten_expressions, x)...)
    else
        return Expression[]
    end
end

function get_coeffs_exponents(f)
    return hcat(collect(Oscar.coefficients(f)), collect(Oscar.exponents(f)))
end


# Translates HC expressions to Oscar polynomials

# ------------------- Input:
# f                   Expression
# varstring           Vector{QQMPolyRingElem}  - Oscar ring gens
# HC_vars             Vector{Variable}         - HC variables
# ------------------- Output:
# f                   QQMPolyRingElem

function HC_to_oscar(f::Expression, varstring::Vector{QQMPolyRingElem}, HC_vars::Vector{Variable})
    E, C = exponents_coefficients(f, HC_vars)
    return sum([C[i] * prod(varstring .^ E[:, i]) for i in 1:length(C)])
end


# ------------------- Input:
# poly                Expression
# vars                Vector{Variable}
# weight              Vector{Int}
# ------------------- Output:
#                     Expression

function get_leading_monomial(poly::Expression, vars::Vector{Variable}, weight::Vector{Int64})

    exponents, coeffs = exponents_coefficients(poly, vars)
    modified_exponents = transpose(exponents) * weight
    leading = minimum(modified_exponents)

    i = findall(x -> x == leading, modified_exponents)[1]
    m = size(exponents)[1]

    return coeffs[i] * prod(j -> vars[j]^exponents[:, i][j], 1:m)
end

# ------------------- Input:
# poly                Expression
# vars                Vector{Variable}
# weight              Vector{Int}
# t                   Variable (Optional)
# ------------------- Output:
#                     Expression

function weight_deformation_for_poly(poly::Expression, vars::Vector{Variable}, weight::Vector{Int64}; t = nothing)

    exponents, coeffs = exponents_coefficients(poly, vars)
    modified_exponents = transpose(exponents) * weight
    leading = minimum(modified_exponents)

    new_exponents = vcat(transpose(modified_exponents .- leading), exponents)
    n = size(new_exponents)[1]
    m = size(new_exponents)[2]

    if t == nothing
        @unique_var t
    end
    new_vars = vcat(t, vars)

    return sum(i -> coeffs[i] * prod(j -> new_vars[j]^new_exponents[:, i][j], 1:n), 1:m)
end


# ------------------- Input:
# parametrization     Vector{Vector{Expression}} or Vector{Expression}  - paramaterization of the map
# random_range        Int                 - random range for valuesfor a point on the variety
# ------------------- Output:
# d                   Int                 - degree of the map

function degree_map(parametrization::Union{Vector{Vector{Expression}}, Vector{Expression}, Vector{Vector{Vector{Expression}}}}; random_range = 100)
    parametrization = normalize_param_groups(parametrization)::Vector{Vector{Expression}}

    l = length(parametrization)

    #homogenize
    @unique_var s[1:l]
    parametrization = [ s[i] .* (parametrization[i]) for i in 1:l]

    #translate to Oscar
    HC_vars = variables(flatten_expressions(parametrization))
    oscar_vars = tuple(HC_vars...)
    R, oscar_vars = polynomial_ring(QQ, string.(HC_vars))
    parametrization_Oscar = vcat([ map(f -> HC_to_oscar(f, oscar_vars, HC_vars), param) for param in parametrization ]...)

    # get a (likely) generic point on the variety
    # We fix the auxiliary homogenizing variables s[i] to 1 to avoid special fibres
    # (e.g. s=0) which can lead to unstable or incorrect degree computations.
    n = length(oscar_vars)
    values = Vector{typeof(QQ(1))}(undef, n)
    s_set = Set(s)
    for (i, v) in enumerate(HC_vars)
        if v in s_set
            values[i] = QQ(1)
        else
            # Avoid zero assignments to reduce the chance of sampling a non-generic
            # point on the image (which can lower the observed fiber cardinality).
            num = rand(vcat(-random_range:-1, 1:random_range))
            den = rand(1:random_range)
            values[i] = QQ(num // den)
        end
    end
    evaluation = [Oscar.evaluate(f, oscar_vars, values) for f in parametrization_Oscar]

    #solve equations for a generic fiber
    I = ideal(parametrization_Oscar - evaluation)
    abs_dec = absolute_primary_decomposition(I)

    return (sum(p[4] for p in abs_dec))
end

# ------------------- Input:
# parametrization     Vector{Vector{Expression}} or Vector{Expression}  - paramaterization of the map
# weight              Vector{Int}
# random_range        Int (Optional)               - random range for valuesfor a point on the variety
# ------------------- Output:
# d                   Int                          - degree of the map

function degree_monomial_map(parametrization::Union{Vector{Vector{Expression}}, Vector{Expression}, Vector{Vector{Vector{Expression}}}}, weight::Vector{Int}; random_range = 100)
    parametrization = normalize_param_groups(parametrization)::Vector{Vector{Expression}}

    l = length(parametrization)

    #cut to monomials
    HC_vars = variables(flatten_expressions(parametrization))
    monom_param = [ [get_leading_monomial(poly, HC_vars, weight) for poly in param] for param in parametrization ]

    return degree_map(monom_param; random_range = random_range)
end


# find number of isolated points in the base locus

# ------------------- Input:
# sagbi               Vector{Expression}  - paramaterization of the map
# ------------------- Output:
# n                   Int                 - number of isolated points in the base locus

function isolated_nsols_baselocus(sagbi::Vector{Expression})

    #translate to Oscar
    HC_vars = variables(sagbi)
    oscar_vars = tuple(HC_vars...)
    R, oscar_vars = polynomial_ring(QQ, string.(HC_vars))
    sagbi_Oscar = map(f -> HC_to_oscar(f, oscar_vars, HC_vars), sagbi)

    I = ideal(sagbi_Oscar)
    abs_dec = absolute_primary_decomposition(I)

    #leave only 0 dimensional components
    abs_dec_points = filter(p -> Oscar.dim(p[3]) == 0, abs_dec)

    if length(abs_dec_points) == 0
        return 0
    else
        return sum(p[4] for p in abs_dec_points)
    end
end

# compute base locus

# ------------------- Input:
# G                   Vector{Expression}
# ------------------- Output:
#

function get_base_locus(G::Vector{Expression})
    return solutions(HomotopyContinuation.solve(G))
end


# compute base locus

# ------------------- Input:
# sys                 System
# ------------------- Output:
# lin_sys
# sagbi
function possible_sagbi(sys::System)

    support, coeffs = support_coefficients(sys)
    vars_sys = variables(sys)

    nonlin_sys_sagbi = map(support, coeffs) do M, c

        unique_elements = unique(c)
        coeffs_monoms = Dict{typeof(c[1]), Expression}()

        for element in unique_elements
            A = M[:, findall(x -> x == element, c)]
            n, m = size(A)
            coeffs_monoms[element] = sum(i -> prod(j -> vars_sys[j]^A[:, i][j], 1:n), 1:m)
        end

        sagbi_part = collect(values(coeffs_monoms))
        fi = transpose(collect(keys(coeffs_monoms))) * sagbi_part

        [fi, sagbi_part]
    end

    n = length(nonlin_sys_sagbi)
    sagbi = union([nonlin_sys_sagbi[i][2] for i in 1:n]...)

    k = length(sagbi)
    @var z[1:k]

    lin_sys = subs([nonlin_sys_sagbi[i][1] for i in 1:n], sagbi => z)

    return lin_sys, sagbi
end


# compute Plücker coordinates of Grassmannian G(k,m)

# ------------------- Input:
# k                   Int
# m                   Int
# ------------------- Output:
# x                   Vector{Variable}     -- HC variables
# sagbi               Vector{Expression}   -- plucker coordinates
# w                   Vector{Int}          -- weight inducing diagonal term order
function get_sagbi_grassmannian(k::Int64, m::Int64)
    @var x[1:(k * (m - k))]

    M = hcat(Matrix(I, k, k), transpose(reshape(x, m - k, k)))

    indx = collect(combinations(1:m, k))
    sagbi = [det(M[:, ind]) for ind in indx]
    w = [i * (m - j + 1) for i in 0:(k - 1) for j in 1:(m - k) ]
    return x, sagbi, w
end
