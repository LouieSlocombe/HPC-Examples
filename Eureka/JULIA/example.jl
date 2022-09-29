
#using Pkg
#Pkg.add("Plots")
using Plots

#using Pkg
#Pkg.add("QuadGK")
using QuadGK


f(x,y,z) = (x^2 + 2y)*z

function integral_calc(y,z)
    arg(x) = f(x,y,z)
    result = quadgk(arg, 1, 10)
    return result
end

result = integral_calc(10,20)
println(result)

x = 1:0.01:10
y = f.(x,10,20)

ENV["GKSwstype"]="100" #set to headless plotting
plot(x, y, label="f(x,10,20)")
plot!(xlab="x", ylab="f(x,10,20)")
savefig("img.png")
