module SimpleED

export fermionicOperators, fermionicOperatorsFunc, fermionicHamiltonianBuilder, measure

function fermionicOperators(sites)
        
    
    f = [[1 0]; [0 -1.]]
    ii = [[1 0]; [0 1.]]
    c = [[0 0]; [1. 0]]
    a = [[0 1.]; [0 0]]
    n = c*a
    
    ci = zeros(Float64, (sites,2^sites,2^sites))
    ai = zeros(Float64, (sites,2^sites,2^sites))
    ni = zeros(Float64, (sites,2^sites,2^sites))
    fi = zeros(Float64, (sites,2^sites,2^sites))
    
    if sites == 1
        ci[1,:,:] = c
        ai[1,:,:] = a
        ni[1,:,:] = n
        fi[1,:,:] = f
    else
        for i in 1:sites
            ci[i,:,:] = kron((f for j in 1:i-1)...,c,(ii for j in i:sites-1)...)
            ai[i,:,:] = kron((f for j in 1:i-1)...,a,(ii for j in i:sites-1)...)
            ni[i,:,:] = kron((ii for j in 1:i-1)...,n,(ii for j in i:sites-1)...)
            fi[i,:,:] = kron((f for j in 1:i)...,(ii for j in i:sites-1)...)
        end
    end
    
    return ci,ai,ni,fi
    end;

function fermionicOperatorsFunc(sites)
    ci,ai,ni,fi = fermionicOperators(sites)
    return (x::Int -> ci[x,:,:], x::Int -> ai[x,:,:], x::Int -> ni[x,:,:], x::Int ->fi[x,:,:])
end

function fermionicHamiltonianBuilder(e,t,U)
    
    
    N = length(e)

    c,a,n,f = fermionicOperatorsFunc(N)
    
    H = zeros(typeof(e[1]*U[1][1]*t[1][1]),(2^N, 2^N))
    
    for i in 1:N
        H .+= e[i]*n(i)
        for j in 1:N
            H .+= t[i,j][1] * c(i) * a(j)
            H .+= U[i,j][1] * n(i) * n(j)
        end
    end
        
    return H
    
end;

function measure(
        observable,
        probabilities,
        ES
    )
    """"""
    
    obsTilde = transpose(ES.vectors)*observable*ES.vectors
    
    output ::typeof(probabilities[1]*observable[1]*ES.vectors[1]) = 0
    
    for (i,P) in enumerate(probabilities)
        output += P * obsTilde[i,i]
    end
    
    return output
    
end;


#function spfermionicOperators(sites)
#        
    
#    f = sparse([[1 0]; [0 -1.]])
#    ii = sparse([[1 0]; [0 1.]])
#    c = sparse([[0 0]; [1. 0]])
#    a = sparse([[0 1.]; [0 0]])
#    n = c*a
    
#    ci = [kron((f for j in 1:i-1)...,c,(ii for j in i:sites-1)...) for i in 1:sites]
#    ai = [kron((f for j in 1:i-1)...,a,(ii for j in i:sites-1)...) for i in 1:sites]
#    ni = [kron((ii for j in 1:i-1)...,n,(ii for j in i:sites-1)...) for i in 1:sites]
#    fi = [kron((f for j in 1:i)...,(ii for j in i:sites-1)...) for i in 1:sites]
    
#    return ci,ai,ni,fi
#    end;

#function spHamiltonianBuilder(e,t,U)
    
#    N = length(e)

#    c,a,n,f = spfermionicOperators(N)
    
#    H = spzeros(Complex{Float64},2^N,2^N)
        
#    for i in 1:N
#        H .+= e[i]*n[i]
#        for j in 1:N
#            H .+= t[i,j][1] * c[i] * a[j]
#            H .+= U[i,j][1] * n[i] * n[j]
#        end
#    end
        
#    return H
    
#end;

end # module
