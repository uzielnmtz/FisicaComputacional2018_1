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


doc"""
newton_numerico(f, x0, h=1e-5, $\epsilon$=1e-5)


Calcula las raices de una funcion mediante el metodo de Newton.
# Ejemplo
julia> f(x) = x^2 - 16


julia> newton_numerico(f, 5)

4.000000000000004
"""
function newton_numerico(f, x0, h=0.001, ϵ=1e-5)
    n = 0
    err = 1
    while ϵ <= err
        n += 1
        x = x0
        x0 = x0 - (f(x0) / der_taylor(f, x0,h))
        err = norm(x - x0)
        if n > 10000
            return x0
        end
    end
    return x0
end


function newton(f, fp, x0, ϵ=1e-5)
    n = 0
    err = 1
    while ϵ <= err
        n += 1
        x = x0
        x0 = x0 - (f(x0) / fp(x0))
        err = norm(x - x0)
        if n > 10000
            return x0
        end
    end
    return x0
end

function newton_sympy(f, x0, ϵ=1e-5)
    using Sympy
    x1 = Sym('x')
    fp = diff(f)
    n = 0
    δ  = 1
    while ϵ <= δ
        n += 1
        x = x0 - (N(f(x0)) / N(fp(x0)))
        δ = norm(x - x0)
        x0 = x
        if n > 10000
            return x0
        end
    end
        return x0
end

function derivada(f, h, x)
    return (f(x + h) - f(x)) / h
end

function der_taylor(f, x, h=1e-5)
    return (f(x + h) - f(x - h)) / (2*h)
end
