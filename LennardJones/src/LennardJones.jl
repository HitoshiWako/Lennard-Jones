module LennardJones
using LinearAlgebra

export Atom,calc_acc,update_x!,update_v!

potential(r,ϵ,σ,p=12,q=6)=4ϵ*((σ/r)^p - (σ/r)^q)
force(r,ϵ,σ,p=12,q=6)=4ϵ/r*(p*(σ/r)^p - q*(σ/r)^q)

mutable struct Atom
    x::Vector{Float64}
    v::Vector{Float64}
    m::Float64
end

function calc_acc(as::Vector{Atom},f::Function,rng=(0.1,1.0))::Vector{Vector{Float64}}
    accs=[]
    for a1 in as
        acc = zeros(size(a1.x))
        for a2 in as
            if a1 != a2 && (rng[1]<abs(norm(a1.x-a2.x))<rng[2])
                acc += f(norm(a1.x-a2.x))/a1.m*(a1.x-a2.x)/norm(a1.x-a2.x)
            end
        end
        push!(accs,acc)
    end
    return accs 
end

function update_x!(as::Vector{Atom},accs::Vector{Vector{Float64}},Δt::Float64)
    for (a, acc) in zip(as,accs)
        a.x +=  a.v*Δt + acc*Δt^2/2
    end
end

function update_v!(as::Vector{Atom},accs::Vector{Vector{Float64}},acc_nexts::Vector{Vector{Float64}},Δt::Float64)
    for (a, acc,acc_next) in zip(as,accs,acc_nexts)
        a.v +=  (acc+acc_next)*Δt/2
    end
end
end # module LennardJones
