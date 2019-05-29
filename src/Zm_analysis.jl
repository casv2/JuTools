module Zm_analysis

using Printf, Plots, JSON, Glob, Statistics
gr(size=(800,500), html_output_format=:png)

export read_Zm, Zm_E, Zm_T, Zm_PT

function read_Zm(filename)
    return JSON.parsefile(filename)
end

function Zm_E(D; return_plot=true)
    p = plot(legend =:false)

    for (index, value) in enumerate(keys(D))
        E_tot, E_kin, E_pot = D[value]["E_tot"], D[value]["E_kin"], D[value]["E_pot"]

        plot!(p,  E_tot)
        plot!(p, E_kin)
        plot!(p, E_pot)
    end

    title!(D["1"]["info"])
    xlabel!("MD step")
    ylabel!("Energy (eV)")

    if return_plot == false
        savefig(@sprintf("E_%s.png", D["1"]["info"]))
    else
        display(p)
    end
end

function Zm_T(D; return_plot=true)
    p = plot(legend =:false)

    for (index, value) in enumerate(keys(D))
        T = D[value]["T"]

        plot!(p, T)
    end

    title!(D["1"]["info"])
    xlabel!("MD step")
    ylabel!("Temperature (K)")

    if return_plot == false
        savefig(@sprintf("T_%s.png", D["1"]["info"]))
    else
        display(p)
    end
end

function Zm_PT(D, delta_step=2, large_data_points=2, eq_length=10; GPa=true, return_plot=true)
    p = scatter(legend =:false)

    colors = ["blue", "red", "green"]

    for (index, value) in enumerate(keys(D))
        T, P = D[value]["T"], D[value]["P"]

        # println(T, P)

        #scatter!(P[eq_length:delta_step:end], T[eq_length:delta_step:end], markersize=1, markerstrokewidth=0, color=colors[index], alpha=0.5)

        large_data_point_step = (length(T)-eq_length)/large_data_points

        mean_T_list, mean_P_list = [], []
        error_T_list, error_P_list = [], []

        for j in 1:large_data_points
            st = round(Int,((j-1)*large_data_point_step + eq_length))
            en = round(Int,(j*large_data_point_step))

            T_slice = Float64.(T[ st : en])
            P_slice = Float64.(P[ st : en])

            mean_T = mean(T_slice)
            mean_P = mean(P_slice)

            error_T = std(T_slice)
            error_P = std(P_slice)

            push!(mean_T_list, mean_T)
            push!(error_T_list, error_T)
            push!(mean_P_list, mean_P)
            push!(error_P_list, error_P)
        end

        if GPa
            mean_P_list = mean_P_list .* 160.21766208
            error_P_list = error_P_list .* 160.21766208
        end

        scatter!(mean_P_list, mean_T_list, xerr=error_P_list, yerr=error_T_list, markersize=5, color=colors[index])
    end

    title!(D["1"]["info"])
    xlabel!("Pressure [GPa]")
    ylabel!("Temperature [K]")

    if return_plot == false
        savefig(@sprintf("PT_%s.png", D["1"]["info"]))
    else
        display(p)
    end
end


end
