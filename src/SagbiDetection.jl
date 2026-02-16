# -------------  Input:
# f              QQMPolyRingElem
# S              QQMPolyRing
# w              Vector{Int64}
# -------------  Output:
# initial_form   QQMPolyRingElem  initial form of a polynomial f w.r.t. w

function initial_form(f::QQMPolyRingElem, S::QQMPolyRing, w::Vector{Int64})
    f_terms = get_coeffs_exponents(f)
    n, _ = size(f_terms)

    if n == 0
        return (f - f)
    end

    Q = MPolyBuildCtx(S)

    term0_coeff, term0_exponent = f_terms[1, 1], f_terms[1, 2]
    push_term!(Q, term0_coeff, term0_exponent)

    d = wDeg(term0_exponent, w)

    for i in 2:n
        term_coeff, term_exponent = f_terms[i, 1], f_terms[i, 2]
        e = wDeg(term_exponent, w)

        if d < e
            finish(Q)
            Q = MPolyBuildCtx(S)
            push_term!(Q, term_coeff, term_exponent)
            d = e
        elseif e == d
            push_term!(Q, term_coeff, term_exponent)
        end
    end
    initial_form = finish(Q)
    return initial_form
end

# ---------------------------------------------- Input:
# G                                              Vector{QQMPolyRingElem}
# R                                              QQMPolyRing
# w                                              Vector{Int64}
# ---------------------------------------------- Output:
# true if G is SAGBI w.r.t. w, false otherwise   Bool

function SagbiCriterion(G::Vector{QQMPolyRingElem}, R::QQMPolyRing, w::Vector{Int})

    n = length(G)
    ord = wdeglex(R, w)

    lt_G = [Oscar.leading_term(G[i], ordering = ord) for i in 1:n]
    #A = [exponent_vector(lt_G[i],1) for i in 1:n]
    #Aw = [wDeg(a, w) for a in A]

    S, vars = Oscar.graded_polynomial_ring(Oscar.QQ, "z" .* string.(1:n))
    g = hom(S, R, lt_G)
    h = hom(S, R, G)
    IA = kernel(g)
    I = kernel(h)

    S1, _ = quo(S, IA)
    S2, _ = quo(S, I)

    return (string(Oscar.hilbert_series(S1)) == string(Oscar.hilbert_series(S2)))
end


# ------------------- Input:
# G                   Vector{QQMPolyRingElem}
# R                   QQMPolyRing
# ------------------- Output:
# list of weights     Vector{Vector{ZZRingElem}}

function weightVectorsRealizingSAGBI(G::Vector{QQMPolyRingElem}, R::QQMPolyRing)
    P = Oscar.newton_polytope(reduce(*, G))
    V = Oscar.vertices(P)

    #normalCones = [Oscar.normal_cone(P, i) for i in 1:length(V)]
    normalCones = maximal_cones(normal_fan(P))

    n = length(V[1])
    id = Oscar.identity_matrix(Oscar.ZZ, n)
    NonnegativeOrth = Oscar.positive_hull(-id)

    for c in normalCones
        intersectedC = Oscar.intersect(c, NonnegativeOrth)
        if Oscar.dim(intersectedC) > 0
            rays = matrix(Oscar.ZZ, Oscar.rays(intersectedC))
            w = ones(Int, size(rays)[1]) * rays

            if (all(wi -> wi < 0, w))
                if (SagbiCriterion(G, R, -Int.(w)))
                    return Int.(w)
                end
            end
        end
    end
    return
end
