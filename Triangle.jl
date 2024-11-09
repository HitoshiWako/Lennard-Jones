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
#---------------------------------------------------------------------------------------------------------------------------
λ = 0.4 # 格子定数
#---------------------------------------------------------------------------------------------------------------------------
a1 = Atom([-λ/√3,λ],[0.0,0.0],m)
a2 = Atom([-λ/√3,-λ],[0.0,0.0],m)
a3 = Atom([λ*2/√3,0.0],[0.0,0.0],m)
as = [a1,a2,a3]

Δt = 1e-6
n = 1e6
samp = 1000

anim = Animation()
@showprogress for i=0:n
    if i%samp == 0
        plt = plot([a1.x[1]],[a1.x[2]],st=:scatter,xlims=(-1.0,1.0),ylims=(-1.0,1.0),label="",aspect_ratio=:equal)
        plot!([a2.x[1]],[a2.x[2]],st=:scatter,label="")
        plot!([a3.x[1]],[a3.x[2]],st=:scatter,label="")
        frame(anim,plt)
    end
    accs = calc_acc(as,f_ar)
    update_x!(as,accs,Δt)
    acc_nexts = calc_acc(as,f_ar)
    update_v!(as,accs,acc_nexts,Δt)

end
gif(anim,"triangle_LJ4.gif",fps=30)