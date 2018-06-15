function read_param(str)
    tmp = match(r"\b[0-9]*.[0-9]*e?[+|-]?[0-9]?\b", str)
    value = parse(Float64, tmp.match)
    return value
end

function input()
    f = open("./d1q3.julia_inp")
    input = Dict{String,Float64}()
    input_list = ["dif", "ar", "br", "cr", "dr", "sigma",
                  "rho0", "rhoinl", "rhoout", "u0", "efield", "nsteps", "diag"]
    for params in input_list
        str = readline(f)
        input[params] = read_param(str) 
    end
end

function init(iseed)
    
end

input()
