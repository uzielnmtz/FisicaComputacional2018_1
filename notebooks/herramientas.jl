function iterar(f, ϵ, x0)
    err = 1
    n = 0
    iters = Float64[]
    numiters= []
    while ϵ <= err
        n += 1
        x = x0
        x0 = f(x0)
        err = norm(x-x0)
        push!(iters, x)
        push!(numiters, n)
    end
    return [numiters iters]
end


function signo(x)
    if x == 0
        return 0
    elseif x > 0
        return 1
    elseif x < 0
        return -1
    end
end


function biseccion(f, a, b)
    iters = []
    if f(a) * f(b) >= 0
        return error
    end
    c = (b + a) / 2
    if f(c) == 0
        return c
    end
    ϵ = 0.000000001
    err = 1
    while ϵ <= err
        if signo(f(a) * f(c)) == -1
            b = c
        else
            a = c
        end
        c = (b + a) / 2
        err = norm(b-a)
        push!(iters, c)
    end
    return iters
end
