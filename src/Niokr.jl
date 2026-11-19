module Niokr
include("ParametersUtils.jl")
include("Components.jl")

using DelimitedFiles: writedlm
using .ParametersUtils: load_from_file, dump_to_file
using .Components: Plus, Minus

export main
function @main(args)
    println("Шаг 1")
    p = open("parameters.toml") do file
        load_from_file(file)
    end
    open("check.toml", "w") do file
        dump_to_file(file, p)
    end

    println("Шаг 2: +")
    A = Plus.A(p, p.V1_0, p.V2_0, 1)
    b = Plus.b(p, p.V1_0, p.V2_0, 1)
    rho = A \ -b
    display_result(A, b, rho)

    println("Шаг 3: -")
    A = Minus.A(p, p.V1_0, p.V2_0, 1)
    b = Minus.b(p, p.V1_0, p.V2_0, 1)
    rho = A \ -b
    display_result(A, b, rho)
end

function display_result(A, b, rho)
    function display_part(name, val)
        println("$name:")
        writedlm(stdout, val)
        println()
    end
    display_part("A", A)
    display_part("b", b)
    display_part("rho", rho)
end

end