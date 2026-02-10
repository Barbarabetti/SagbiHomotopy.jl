using SagbiHomotopy
using HomotopyContinuation

function oscillators(N, M)
    @var u[1:M, 1:N], v[1:M,1:N]
    @var z[1:2*M*N, 1:M+3+N-1]
    
    polynomials = Matrix{Expression}(undef, 2*M*N, M+3+N-1);
    
    for k in 1:N   
        for i in 1:M

            polynomials[2*i - 1 + 2*M*(k-1),1] = (u[1,1])^0
            polynomials[2*i - 1 + 2*M*(k-1),2] = u[i, k]
            polynomials[2*i - 1 + 2*M*(k-1),3] = v[i, k]

            polynomials[2*i + 2*M*(k-1),1] = (u[1,1])^0
            polynomials[2*i + 2*M*(k-1),2] = u[i, k]
            polynomials[2*i + 2*M*(k-1),3] = v[i, k]

            for j in 1:M
                polynomials[2*i - 1 + 2*M*(k-1),j+3]= u[i, k]*(u[j, k]^2 + v[j, k]^2)
            end

            for j in 1:(N-1)
                lst = filter(j -> j != k, collect(1:N))
                polynomials[2*i - 1 + 2*M*(k-1),j+M+3]= v[i, lst[j]]
            end

            for j in 1:M
                polynomials[2*i + 2*M*(k-1),j+3]= v[i, k]*(u[j, k]^2 + v[j, k]^2)
            end

            for j in 1:(N-1)
                lst = filter(j -> j != k, collect(1:N))
                polynomials[2*i + 2*M*(k-1),j+M+3]= u[i, lst[j]]
            end
        end
    end

    sagbi = [polynomials[i, :] for i in 1:2*M*N]
    Sagbi = []
    for i = 1:2:(length(sagbi)-1)
        s = unique(vcat(sagbi[i],sagbi[i+1]))
        push!(Sagbi,s)
    end

    c = rand(Float64, 2*M*N, M+3+N-1);
    lin_sys = [[sum(z[i,:].* c[i,:])] for i in 1:2*M*N];

    Lin_sys = []
    for i = 1:2:(length(lin_sys)-1)
        l = unique(vcat(lin_sys[i], subs(lin_sys[i+1], z[i+1,1]=>z[i,1], z[i+1,2]=>z[i,2], z[i+1,3]=> z[i,3] )))
        push!(Lin_sys,l)
    end
    
    return(Lin_sys, Sagbi)
end

#run the function for the first time without storing the time
F1, sagbi1 = oscillators(1, 5);
w1 = vcat(fill(-2, 5), fill(-1, 5));

sagbi_homotopy(F1, sagbi1, weight = w1, degreeCheck = false);

S = vcat([subs(F1[i], variables(F1[i])=>sagbi1[i]) for i in 1:length(F1)]...)
solve(S)
N = 5
M = 5


for s=1:N
    for r=1:M
        F, sagbi = oscillators(s, r)
        w = vcat(fill(-2, s*r), fill(-1, s*r))

        result, elapsed_time, allocations, _, _ = @timed sagbi_homotopy(F, sagbi, weight = w, degreeCheck = false, varyLinearPart = true)
        #open("Oscillators/oscillators_timing_log_sagbi.txt", "a") do io
        #    println(io, "N = ", s, ", ", "M = ", r)
        #    println(io, "degree: ", length(result[2]))
        #    println(io, "Execution time: ", elapsed_time, " seconds")
        #    println(io, "Memory allocated: ", allocations, " bytes")
        #end 
    end
end