using Pkg
Pkg.activate("LennardJones")

using Plots
using Genie,Genie.Renderer.Json,Genie.Requests
using Format

using LennardJones
#----------------- utility function ----------------------
function seriesfilename(prefix::String,suffix::String,dir::String=".")::String
    filenames = readdir(dir)
    i = 1
    while true
        fn = format("{:s}{:06d}{:s}",prefix,i,suffix) 
        if fn ∉ filenames
            touch(dir*"/"*fn)
            return fn
        end
        i += 1
    end
end
#----------------- Parameters ---------------------------
Δt = 1e-6 #　10^-12sec(ps)
n = 1e6
samp = 10000
#---------------- Simulation Function -------------------
function verlet(fname::String,as::Vector{Atom},func::Function,
                    xlims::Tuple{Float64,Float64},ylims::Tuple{Float64,Float64})
    
    anim = Animation()
    for i=0:n
        if i%samp == 0
            xs = [a.x[1] for a in as]
            ys = [a.x[2] for a in as]
            time = i/1000
            plt = plot(xs,ys,st=:scatter,xlims=xlims,ylims=ylims,label="",
                       aspect_ratio=:equal,title="$time[ns]")
            frame(anim,plt)
        end
        accs = calc_acc(as,func)
        update_x!(as,accs,Δt)
        acc_nexts = calc_acc(as,func)
        update_v!(as,accs,acc_nexts,Δt)
    end
    gif(anim,"result/"*fname,fps=30)
end

function singleparticle(fname::String,r::Float64,ϵ::Float64,σ::Float64,m::Float64)
    rs = 0.2:0.01:1.0
    a = Atom([r],[0.0],m)
    anim = Animation()
    for i=1:n
        if i%samp == 0
            plt = plot(rs,LennardJones.potential.(rs,ϵ,σ),xlims=(0.2,1.0),ylims=(-2,2),label="",title="Lennard Jones Potential",xlabel="distance[nm]",ylabel="Energy[zJ]")
            plot!([a.x[1]],[LennardJones.potential(a.x[1],ϵ,σ)],st=:scatter,label="")
            frame(anim,plt)
        end
        acc = LennardJones.force(a.x[1],ϵ,σ)/a.m
        a.x[1] += a.v[1]*Δt+ acc*Δt^2/2
        acc_next = LennardJones.force(a.x[1],ϵ,σ)/a.m
        a.v[1] += (acc+acc_next)/2*Δt
    end
    gif(anim,"result/"*fname,fps=30)
end

#----------------- Routing ------------------------------
route("/singleparticle",method=POST) do
    msg = jsonpayload()
    r = msg["r0"]
    ϵ = msg["eps"]
    σ = msg["sigma"]
    m = msg["m"]

    fname = seriesfilename("SingleParticle",".gif","result")
    Threads.@spawn singleparticle(fname,r,ϵ,σ,m)
    (:filename => fname) |> json
end

route("/doubleparticles",method=POST) do
    msg = jsonpayload()
    r = msg["r0"]
    ϵ = msg["eps"]
    σ = msg["sigma"]
    m = msg["m"]
    
    fname = seriesfilename("DoubleParticles",".gif","result")
    as = vec([Atom([-r/2,0],[0.0,0.0],m),Atom([ r/2,0],[0.0,0.0],m)])
    func(r) = LennardJones.force(r,ϵ,σ)
    Threads.@spawn verlet(fname,as,func,(-r,r),(-r,r))
    (:filename => fname) |> json
end

route("/triangle",method=POST) do
    msg = jsonpayload()
    λ = msg["lam"]
    ϵ = msg["eps"]
    σ = msg["sigma"]
    m = msg["m"]
    
    fname = seriesfilename("Triangle",".gif","result")
    as = vec([Atom([-λ/2,-λ*√3/4],[0.0,0.0],m),
              Atom([ λ/2,-λ*√3/4],[0.0,0.0],m),
              Atom([   0, λ*√3/4],[0.0,0.0],m)])
    func(r) = LennardJones.force(r,ϵ,σ)
    Threads.@spawn verlet(fname,as,func,(-λ,λ),(-λ,λ))
    (:filename => fname) |> json
end
route("/lattice",method=POST) do
    msg = jsonpayload()
    λ = msg["lam"]
    cycle = msg["cycle"]
    ϵ = msg["eps"]
    σ = msg["sigma"]
    m = msg["m"]
    
    fname = seriesfilename("Lattice",".gif","result")
    as = vec([Atom([x,y],[0.0,0.0],m) for x=-0.0:λ:λ*(cycle-1),y=0.0:λ:λ*(cycle-1)])
    func(r) = LennardJones.force(r,ϵ,σ)
    Threads.@spawn verlet(fname,as,func,(-λ,λ*cycle),(-λ,λ*cycle))
    (:filename => fname) |> json
end

route("/trianglelattice",method=POST) do
    msg = jsonpayload()
    λ = msg["lam"]
    row = msg["row"]
    col = msg["col"]
    ϵ = msg["eps"]
    σ = msg["sigma"]
    m = msg["m"]
    
    fname = seriesfilename("TriangleLattice",".gif","result")
    as = vec([Atom([x,y],[0.0,0.0],m) for x=0.0:λ:λ*(col-1),y=0.0   :λ*√3:λ*√3*(row-1)/2])
    bs = vec([Atom([x,y],[0.0,0.0],m) for x=λ/2:λ:λ*(col-1),y=λ*√3/2:λ*√3:λ*√3*(row-1)/2])
    append!(as,bs)
    func(r) = LennardJones.force(r,ϵ,σ)
    Threads.@spawn verlet(fname,as,func,(-λ,λ*col),(-λ,λ*(√3*row/2)))
    (:filename => fname) |> json
end
# ---------------- Run Server -----------------------------
up(async=false)
# ---------------- クライアント側で以下を実行 ------------------
# curl -X POST -H "Content-Type: application/json" \
#      -d '{"r0":0.8, "eps":1.67, "sigma":0.34, "m":6.63e-3 }' \
#      "http://127.0.0.1:8000/singleparticle"