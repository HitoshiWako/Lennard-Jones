#---------------------------------------------------------------------------------------------------------------------------
using Pkg
Pkg.activate("LennardJones")
Pkg.add("Plots")
Pkg.add("ProgressMeter")

using LennardJones

using Plots
using ProgressMeter
#---------------------------------------------------------------------------------------------------------------------------
# Parameter
ϵ = 1.67
σ = 0.34
m = 6.63e-3

f_ar(r) = LennardJones.force(r,ϵ,σ)

a = Atom([0.8],[0.0],m)
Δt = 1e-6
n = 1e6
samp = 1000

rs = 0.2:0.01:1.0

anim = Animation()
@showprogress for i=1:n
    if i%samp == 0
        plt = plot(rs,LennardJones.potential.(rs,ϵ,σ),xlims=(0.2,1.0),ylims=(-2,2),label="",title="Lennard Jones Potential",xlabel="distance[nm]",ylabel="Energy[zJ]")
        plot!([a.x[1]],[LennardJones.potential(a.x[1],ϵ,σ)],st=:scatter,label="")
        frame(anim,plt)
    end
    acc = f_ar(a.x[1])/a.m
    a.x[1] += a.v[1]*Δt+ acc*Δt^2/2
    acc_next = f_ar(a.x[1])/a.m
    a.v[1] += (acc+acc_next)/2*Δt
end
gif(anim,"result/SingleParticle.gif",fps=30)