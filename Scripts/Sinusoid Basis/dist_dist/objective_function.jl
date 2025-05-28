using Ipopt
using LinearAlgebra
using StaticArrays

function cost_function(model, slack)
    @objective(model, Min, slack)
end