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
row = 10::Int64
col = 10::Int64
#---------------------------------------------------------------------------------------------------------------------------
as = vec([Atom([(x-1)*λ,(y-1)*√3λ],[0.0,0.0],m) for x=1:col, y=1:div(row,2)+mod(row,2)])
bs = vec([Atom([(x-1/2)*λ,(y-1/2)*√3λ],[0.0,0.0],m) for x=1:col-1, y=1:div(row,2)])
append!(as,bs)
Δt = 1e-6
n = 1e6
samp = 1000

anim = Animation()
@showprogress for i=0:n
    if i%samp == 0
        xs = [a.x[1] for a in as]
        ys = [a.x[2] for a in as]
        time = i/samp
        plt = plot(xs,ys,st=:scatter,xlims=(-λ/2,(col-1/2)*λ),ylims=(-λ/2,(row-1/2)*√3/2*λ),label="",aspect_ratio=:equal,title="$time[ns]")
        frame(anim,plt)
    end
    accs = calc_acc(as,f_ar)
    update_x!(as,accs,Δt)
    acc_nexts = calc_acc(as,f_ar)
    update_v!(as,accs,acc_nexts,Δt)
end
gif(anim,"triangle_lattice3.gif",fps=30)
#---------------------------------------------------------------------------------------------------------------------------