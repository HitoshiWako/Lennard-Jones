# REPLにて以下を実行
# using HTTP
# HTTP.request(
#     "POST",
#     "http://127.0.0.1:8000/solve",
#     [("Content-Type", "application/json")],
#     """{"r0":0.8,"eps":1.67,"sigma":0.34,"m":6.63e-3}"""
# )
#---------------------------------------------------------------------------------------------------------------------------
using Pkg
Pkg.activate("LennardJones")
#Pkg.add("Plots")
#Pkg.add("ProgressMeter")
#Pkg.add("Format")
#Pkg.add("Genie")

using Plots
#using ProgressMeter
using Genie,Genie.Renderer.Json,Genie.Requests
using Format

using LennardJones
#---------------------------------------------------------------------------------------------------------------------------
# Parameter
#ϵ = 1.67
#σ = 0.34
#m = 6.63e-3

#f_ar(r) = LennardJones.force(r,ϵ,σ)

# utility function
function seriesfilename(prefix,suffix,dir=".")
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


Δt = 1e-6
n = 1e6
samp = 1000

function singleparticle(fname,r,ϵ,σ,m)
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

route("/singleparticle",method=POST) do
    msg = jsonpayload()
    fname = seriesfilename("SingleParticle",".gif","result")
    r = msg["r0"]
    ϵ = msg["eps"]
    σ = msg["sigma"]
    m = msg["m"]
    Threads.@spawn singleparticle(fname,r,ϵ,σ,m)
    (:filename => fname) |> json
end
up(async = false)
