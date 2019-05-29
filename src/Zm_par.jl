module Zm_par

using Printf, JuLIP, Distributions, LinearAlgebra, ASE

export Stationary, VelocityVerlet, Zmethod, MaxwellBoltzmann_scale, MaxwellBoltzmann, kB, fs

kB = 8.617330337217213e-05 #units taken from Python ASE
fs = 0.09822694788464063

function MaxwellBoltzmann(at, temp)
    d = Normal()
    M = reshape(collect(rand(d, 3*length(at.M))), (length(at.M),3))
    M2 = M .* sqrt.(at.M .* temp)
    set_momenta!(at, collect(transpose(M2)))

    return at
end

function Stationary(at)
    p0 = sum(at.P)
    mtot = sum(at.M)
    v0 = p0 / mtot

    momenta = [at.M[i] * v0 for i in 1:length(at.M)]
    set_momenta!(at, collect(at.P - momenta))

    return at
end

function VelocityVerlet(IP, at, dt)
    V = at.P ./ at.M
    A = forces(IP, at) ./ at.M

    set_positions!(at, at.X + (V .* dt) + (.5 * A * dt^2))

    nA = forces(IP, at) ./ at.M
    nV = V + (.5 * (A + nA) * dt)
    set_momenta!(at, collect(nV .* at.M))

    return at
end

function MaxwellBoltzmann_scale(at, temp)
    d = Normal()
    M = reshape(collect(rand(d, 3*length(at.M))), (length(at.M),3))

    N = []

    for j in 1:1000
        d = Normal()
        M = reshape(collect(rand(d, 3*length(at.M))), (length(at.M),3))
        n = 0
        for i in 1:1:length(M[:,1])
            n += norm(M[i,:])
        end
        push!(N,n)
    end

    d = Normal()
    M = reshape(collect(rand(d, 3*length(at.M))), (length(at.M),3))
    Mnm = mean(N)

    for i in 1:10000
        n = 0
        for i in 1:length(M[:,1])
            n += norm(M[i,:])
        end
        if n < Mnm
            M = M*1.0001
        else
            M = M*0.9999
        end
    end

    M2 = M .* sqrt.(at.M .* temp)
    set_momenta!(at, transpose(M2))

    return at #, Ml, Mnm, Sl, M
end

function Zmethod(IP, at, nsteps, dt, A, N, save_config, element)
    E0 = energy(IP, at)

    m = at.M

    E_tot = zeros(nsteps)
    E_pot = zeros(nsteps)
    E_kin = zeros(nsteps)
    P = zeros(nsteps)
    T = zeros(nsteps)

    al = []

    for i in 1:nsteps

        at = VelocityVerlet(IP, at, dt * fs)
        Ek = ((0.5 * sum(at.M) * norm(at.P ./ at.M)^2)/length(at.M)) / length(at.M)
        Ep = (energy(IP, at) - E0) / length(at.M)
        E_tot[i] = Ek + Ep
        E_pot[i] = Ep
        E_kin[i] = Ek
        T[i] = Ek / (1.5 * kB)
        P[i] = -tr(stress(IP, at))/3.0

        v = at.P ./ m
        C = A/norm(v)

        set_momenta!(at, collect((v + C*v) .* m))

        if i % save_config == 0
            push!(al, deepcopy(at))
        end
    end

    return E_tot, E_pot, E_kin, P, T, al
end

end
