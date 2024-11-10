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
a1 = Atom([0.35],[0.0],m)
a2 = Atom([-0.35],[0.0],m)
as = [a1,a2]
Δt = 1e-6
n = 1e6
samp = 1000

anim = Animation()
@showprogress for i=0:n
    if i%samp == 0
        plt = plot([as[1].x],[0.0],st=:scatter,xlims=(-0.5,0.5),ylims=(-0.5,0.5),label="")
        plot!([as[2].x],[0.0],st=:scatter,label="")
        frame(anim,plt)
    end

    accs = calc_acc(as,f_ar)
    update_x!(as,accs,Δt)
    acc_nexts = calc_acc(as,f_ar)
    update_v!(as,accs,acc_nexts,Δt)

#    acc1 = f_ar(a1.x-a2.x)/a1.m
#    acc2 = f_ar(a2.x-a1.x)/a2.m
#    a1.x += a1.v*Δt+acc1*Δt^2/2
#    a2.x += a2.v*Δt+acc2*Δt^2/2
#    acc_next1 = f_ar(a1.x-a2.x)/a1.m
#    acc_next2 = f_ar(a2.x-a1.x)/a2.m
#    a1.v += (acc1+acc_next1)/2*Δt
#    a2.v += (acc2+acc_next2)/2*Δt

end
gif(anim,"double_LJ2.gif",fps=30)