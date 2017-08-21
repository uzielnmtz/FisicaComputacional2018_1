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
