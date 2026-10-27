# =============================================================================
# safespace.jl  —  reusable Simcyp-style safe-space figure builder.
#
# `include`d (not a module) AFTER `using .BEUtils`, CSV, DataFrames, Statistics,
# CairoMakie. Used by 03_safe_space_plot.jl (ER scenario) and
# 03b_safe_space_plot_IR.jl (IR high-variability scenario).
#
# build_safespace(pk, out; cv_map, tag, scenario, ...) writes one figure per
# endpoint (Cmax, AUC_inf, AUC_0-t, joint), each with n=48 & n=86 panels, and
# returns a DataFrame of safe-space boundaries (power >= power_target).
# =============================================================================

"""
    build_safespace(pk, out; cv_map, tag, scenario, n_trial, n_subj_list,
                    power_target = 0.80, dx = 0.05) -> DataFrame

`cv_map`  : Dict endpoint => within-subject CV (keys "Cmax","AUC_inf","AUC_tEnd").
`tag`     : filename token ("" -> fig_safespace_<ep>.png; "IR" -> ..._IR_<ep>.png).
`scenario`: label shown in the figure super-title (e.g. "ER CV", "IR CV").
"""
function build_safespace(pk, out; cv_map, tag = "", scenario = "",
                         n_trial = N_TRIAL_DEFAULT, n_subj_list = N_SUBJ_LIST,
                         power_target = 0.80, dx = 0.05)
    anchors = sort(unique(pk.T50))
    ids     = sort(unique(pk.IndividualId))
    grid    = collect(range(anchors[1], anchors[end]; step = dx))
    ref_j   = findfirst(==(REF_T50), anchors)
    pref    = isempty(tag) ? "" : tag * "_"
    sclab   = isempty(scenario) ? "" : "  ($scenario)"

    metric_matrix(metric) = begin
        M = Matrix{Union{Missing,Float64}}(missing, length(ids), length(anchors))
        for (j, t) in enumerate(anchors)
            d = sort(filter(:T50 => ==(t), pk), :IndividualId)
            M[:, j] = d[!, metric]
        end
        M
    end
    interp_log(vals, t) = begin
        ok = .!ismissing.(vals); sum(ok) < 2 && return missing
        x = anchors[ok]; v = Float64.(vals[ok])
        t <= x[1]   && return v[1]
        t >= x[end] && return v[end]
        k = searchsortedlast(x, t); (t == x[k]) && return v[k]
        w = (log(t) - log(x[k])) / (log(x[k+1]) - log(x[k]))
        exp((1 - w) * log(v[k]) + w * log(v[k+1]))
    end
    curve_single(M, cv, nss) = begin
        rv = M[:, ref_j]
        [vbe_replicate([interp_log(M[i, :], t) for i in eachindex(ids)], rv;
                       n_subj = nss, n_trial = n_trial, wsv_cv = cv).pass_rate for t in grid]
    end
    anchor_single(M, cv, nss) = begin
        rv = M[:, ref_j]
        [vbe_replicate(M[:, j], rv; n_subj = nss, n_trial = n_trial, wsv_cv = cv).pass_rate
         for j in eachindex(anchors)]
    end
    curve_joint(MC, MA, nss) = begin
        rC = MC[:, ref_j]; rA = MA[:, ref_j]
        [vbe_replicate_joint([interp_log(MC[i, :], t) for i in eachindex(ids)], rC,
                             [interp_log(MA[i, :], t) for i in eachindex(ids)], rA;
                             n_subj = nss, n_trial = n_trial,
                             cv_C = cv_map["Cmax"], cv_A = cv_map["AUC_tEnd"]).pass_rate
         for t in grid]
    end
    anchor_joint(MC, MA, nss) = begin
        rC = MC[:, ref_j]; rA = MA[:, ref_j]
        [vbe_replicate_joint(MC[:, j], rC, MA[:, j], rA; n_subj = nss, n_trial = n_trial,
                             cv_C = cv_map["Cmax"], cv_A = cv_map["AUC_tEnd"]).pass_rate
         for j in eachindex(anchors)]
    end
    safe_bounds(power) = begin
        above = power .>= power_target
        ri = argmin(abs.(grid .- REF_T50)); lo = ri; hi = ri
        while lo > 1 && above[lo-1]; lo -= 1; end
        while hi < length(grid) && above[hi+1]; hi += 1; end
        (grid[lo], grid[hi])
    end
    draw_panel!(ax, power, anchor_pow; show_ylabel = true) = begin
        xlims!(ax, anchors[1], anchors[end]); ylims!(ax, -0.02, 1.02)
        for (k, t) in enumerate(grid)
            p = power[k]
            col = p >= power_target ?
                  (:seagreen,  0.15 + 0.45 * (p - power_target) / (1 - power_target)) :
                  (:firebrick, 0.10 + 0.18 * (power_target - p) / power_target)
            poly!(ax, Point2f[(t-dx/2,-0.02),(t+dx/2,-0.02),(t+dx/2,1.02),(t-dx/2,1.02)],
                  color = col, strokewidth = 0)
        end
        lo, hi = safe_bounds(power)
        hlines!(ax, [0.05, power_target]; linestyle = :dot, color = :gray30)
        vlines!(ax, [REF_T50]; linestyle = :dash, color = :gray40)
        vlines!(ax, [lo, hi]; linestyle = :dot, color = :black)
        lines!(ax, grid, power; color = :black, linewidth = 2.5)
        scatter!(ax, anchors, anchor_pow; color = :white, strokecolor = :black,
                 strokewidth = 1.5, markersize = 12)
        text!(ax, (lo+hi)/2, 0.45;
              text = "SAFE SPACE\n$(round(lo,digits=2))–$(round(hi,digits=2)) h",
              align = (:center,:center), fontsize = 13, color = :darkgreen, font = :bold)
        show_ylabel || hideydecorations!(ax, grid = false)
        (lo, hi)
    end
    endpoint_figure(label, token, curve_fun, anchor_fun, bnd) = begin
        fig = Figure(size = (1150, 500), fontsize = 16)
        Label(fig[0, 1:2], "Metformin XR safe space — $label$sclab  (vs T50=3.0 h ref, $n_trial trials)";
              fontsize = 18, font = :bold)
        for (col, nss) in enumerate(n_subj_list)
            ax = Axis(fig[1, col], xlabel = "Test formulation T50 (h)",
                      ylabel = "P(conclude BE)  —  power", title = "n = $nss subjects",
                      xlabelsize = 17, ylabelsize = 17, titlesize = 17,
                      xticklabelsize = 14, yticklabelsize = 14,
                      xticks = anchors)   # a tick at every simulated T50
            lo, hi = draw_panel!(ax, curve_fun(nss), anchor_fun(nss); show_ylabel = col == 1)
            push!(bnd, (scenario = isempty(scenario) ? "ER" : scenario,
                        endpoint = label, n_subj = nss, ss_lo = lo, ss_hi = hi))
        end
        save(joinpath(out, "fig_safespace_$(pref)$(token).png"), fig, px_per_unit = 2)
    end

    Cmat = metric_matrix("Cmax"); Imat = metric_matrix("AUC_inf"); Amat = metric_matrix("AUC_tEnd")
    bnd  = DataFrame()
    endpoint_figure("Cmax", "Cmax",
        nss -> curve_single(Cmat, cv_map["Cmax"], nss),
        nss -> anchor_single(Cmat, cv_map["Cmax"], nss), bnd)
    endpoint_figure("AUC_inf", "AUCinf",
        nss -> curve_single(Imat, cv_map["AUC_inf"], nss),
        nss -> anchor_single(Imat, cv_map["AUC_inf"], nss), bnd)
    endpoint_figure("AUC_0-t", "AUC0t",
        nss -> curve_single(Amat, cv_map["AUC_tEnd"], nss),
        nss -> anchor_single(Amat, cv_map["AUC_tEnd"], nss), bnd)
    endpoint_figure("Joint (Cmax & AUC_0-t)", "joint",
        nss -> curve_joint(Cmat, Amat, nss),
        nss -> anchor_joint(Cmat, Amat, nss), bnd)
    bnd
end
