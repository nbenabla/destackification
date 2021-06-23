module BerghAlgorithms
include("StackyFan.jl")
include("AlgC.jl")
#Right now, StackyFan.jl contains the module's code -- should bring StackyFan here
export BerghA

function __init__()

    println("         ")
    println("   ^     ")
    println("  / \\    | WARNING! The BerghAlgorithms module is new and key features are still being implemented and tested.")
    println(" / ! \\   | As such, please be patient as bugs and mathematical inaccuracies are resolved.")
    println("-------  | If you find an issue, please consider reporting it to htalbott@umich.edu.")

end

# Write your package code here.


end
