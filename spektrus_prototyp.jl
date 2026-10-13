### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 81eb466a-8faa-420b-a51e-36781af7868b
md"""
# Spektruś.jl — prototyp

Eksplorator widm emisyjnych na danych z NIST ASD (`spectrum.csv`). Wybierz pierwiastek i zakres,
przełącz próżnię/powietrze i skalę natężeń. Pseudokolory: UV = fiolet, podczerwień = ciemna czerwień.
"""

# ╔═╡ 9735532e-8621-4f11-bec1-f557068e179e
begin
    import Pkg
    Pkg.activate(@__DIR__)   # środowisko projektu (Project.toml obok notebooka)
    using CairoMakie, CSV, DataFrames, PlutoUI
end

# ╔═╡ b3c367fa-571f-4406-b32b-e3bfc541f3f0
begin
    "Długość fali z próżni na powietrze (wzór Edlena/Mortona), tylko dla 200–2000 nm."
    function vac2air(λ)
        (200 <= λ <= 2000) || return λ
        s2 = (1000 / λ)^2
        n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
        return λ / n
    end

    "Przybliżony kolor światła widzialnego (algorytm Brutona); pseudokolory dla UV i IR."
    function wl2rgb(λ)
        if λ < 380
            return (0.62, 0.45, 0.95)
        elseif λ > 780
            return (0.80, 0.25, 0.20)
        end
        r, g, b = if λ < 440
            (-(λ - 440) / 60, 0.0, 1.0)
        elseif λ < 490
            (0.0, (λ - 440) / 50, 1.0)
        elseif λ < 510
            (0.0, 1.0, -(λ - 510) / 20)
        elseif λ < 580
            ((λ - 510) / 70, 1.0, 0.0)
        elseif λ < 645
            (1.0, -(λ - 645) / 65, 0.0)
        else
            (1.0, 0.0, 0.0)
        end
        f = λ < 420 ? 0.3 + 0.7 * (λ - 380) / 40 :
            λ > 700 ? 0.3 + 0.7 * (780 - λ) / 80 : 1.0
        γ = 0.8
        return ((f * r)^γ, (f * g)^γ, (f * b)^γ)
    end
end

# ╔═╡ 9b1b3791-b66e-4828-9a43-f4d2b9aa7ae7
begin
    str_cols = Dict(c => String for c in
        [:symbol, :Aki, :conf_i, :term_i, :J_i, :conf_k, :term_k, :J_k])
    df = CSV.read(joinpath(@__DIR__, "spectrum.csv"), DataFrame; types = str_cols)
    elements = sort(combine(groupby(df, [:Z, :symbol]), nrow => :n), :Z)
    element_options = [r.Z => "$(r.symbol) ($(r.Z))" for r in eachrow(elements)]
    have_model = all(c -> c in propertynames(df), [:Aki_s, :Ek_cm1, :gk])
    warn_txt = have_model ? "" :
        "⚠️ Brak kolumn do modelu LTE: uruchom ponownie fetch_nist_lines.jl i prepare_lines.jl."
    md"Wczytano **$(nrow(df))** linii dla **$(nrow(elements))** pierwiastków. $(warn_txt)"
end

# ╔═╡ 7b9c7753-081d-4d70-ab3c-f20def9c4a5d
md"""
**Pierwiastek:** $(@bind elem PlutoUI.Select(element_options; default = 11))

**Od (nm):** $(@bind lo PlutoUI.Slider(100:1:2500; default = 380, show_value = true))
**Do (nm):** $(@bind hi PlutoUI.Slider(100:1:2500; default = 780, show_value = true))

**Długość fali:** $(@bind medium PlutoUI.Select(["air" => "powietrze", "vac" => "próżnia"]))
   **Skala natężeń:** $(@bind scale_mode PlutoUI.Select(["log" => "logarytmiczna", "lin" => "liniowa"]))
   $(@bind ion2 PlutoUI.CheckBox(default = false)) pokaż także jony II

**Rozdzielczość** [nm] (0 = wyłączona; linie w jednym przedziale są łączone, natężenia sumowane): $(@bind res PlutoUI.Slider(0:0.1:5; default = 0.0, show_value = true))

**Próg widoczności** (ile rzędów wielkości poniżej najsilniejszej linii): $(@bind decades PlutoUI.Slider(0.5:0.1:6; default = 4.0, show_value = true))

**Źródło natężeń:** $(@bind intens_src PlutoUI.Select(["lte" => "model LTE: g·A·ν·exp(−E/kT)", "nist" => "NIST (jednostki umowne)"]))
   **T [K]:** $(@bind temp PlutoUI.Slider(1000:100:10000; default = 6000, show_value = true))

$(@bind show_exp PlutoUI.CheckBox(default = false)) zaznacz linie ze sprawozdania (spektroskop pryzmatyczny, tylko Na)
"""

# ╔═╡ 79cead40-7703-4893-937f-a1b0b4e0b070
begin
    λlo = lo
    λhi = max(hi, lo + 1)
    sub = df[(df.Z .== elem) .& (df.ion .<= (ion2 ? 2 : 1)), :]
    sub = transform(sub, :wl_vac_nm => ByRow(x -> medium == "air" ? vac2air(x) : x) => :λ)
    if intens_src == "lte" && have_model
        # emisja optycznie cienka w LTE: g_k * A_ki * nu * exp(-E_k / kT), c2 = hc/k = 1.438777 cm*K
        sub = dropmissing(sub, [:Aki_s, :Ek_cm1, :gk])
        sub.I = [r.gk * r.Aki_s * (1e7 / r.wl_vac_nm) * exp(-1.438777 * r.Ek_cm1 / temp)
                 for r in eachrow(sub)]
        sub = sub[sub.I .> 0, :]
    else
        sub.I = sub.intensity
    end
    if res > 0 && nrow(sub) > 0
        # jak pryzmat/piksel: linie w jednym przedziale o szerokości res zlewają się w jedną,
        # natężenia sumujemy, a położenie to średnia ważona natężeniem
        sub.bin = floor.(Int, sub.λ ./ res)
        sub = combine(groupby(sub, :bin)) do g
            w = g.I
            k = argmax(w)
            (symbol = g.symbol[k], ion = g.ion[k],
             λ = sum(g.λ .* w) / sum(w), wl_vac_nm = sum(g.wl_vac_nm .* w) / sum(w),
             I = sum(w), n_linii = nrow(g), intensity = maximum(g.intensity),
             Aki = g.Aki[k], conf_i = g.conf_i[k], term_i = g.term_i[k], J_i = g.J_i[k],
             conf_k = g.conf_k[k], term_k = g.term_k[k], J_k = g.J_k[k])
        end
    else
        sub.n_linii = ones(Int, nrow(sub))
    end
    inwin = sub[(sub.λ .>= λlo) .& (sub.λ .<= λhi), :]
    lmin, lmax = isempty(sub) ? (0.0, 1.0) : extrema(log10.(sub.I))
    imax = isempty(sub) ? 1.0 : maximum(sub.I)
    cutoff = lmax - decades                       # próg w skali log10
    vis = inwin[log10.(inwin.I) .>= cutoff, :]
    relheight(I) = scale_mode == "log" ?
        clamp((log10(I) - cutoff) / decades, 0.04, 1.0) : I / imax
    nothing
end

# ╔═╡ fc11e62f-5c5a-40cc-ae64-71d9f7117dbd
begin
    src_label = intens_src == "lte" ? "model LTE (T = $(temp) K)" : "NIST"
    Markdown.parse("Widocznych linii: **$(nrow(vis))** z $(nrow(inwin)) w zakresie " *
        "(łącznie dla pierwiastka: $(nrow(sub))). Źródło natężeń: **$(src_label)**; " *
        "zakres $(round(10^lmin; sigdigits = 2)) – $(round(10^lmax; sigdigits = 2)) (jednostki względne).")
end

# ╔═╡ 533e7dde-9e80-4606-abc8-5747a5d7d29c
begin
    fig = Figure(size = (1000, 480), backgroundcolor = RGBf(0.09, 0.09, 0.11))
    axstyle = (backgroundcolor = :black,
               xlabelcolor = :white, ylabelcolor = :white,
               xticklabelcolor = :white, yticklabelcolor = :white,
               xtickcolor = :white, ytickcolor = :white,
               leftspinecolor = :gray40, bottomspinecolor = :gray40,
               topspinevisible = false, rightspinevisible = false,
               xgridvisible = false, ygridvisible = false)
    ax1 = Axis(fig[1, 1]; height = 90,
               yticksvisible = false, yticklabelsvisible = false,
               xticksvisible = false, xticklabelsvisible = false, axstyle...)
    ax2 = Axis(fig[2, 1];
               xlabel = "λ [nm, $(medium == "air" ? "powietrze" : "próżnia")]",
               ylabel = "natężenie względne", axstyle...)

    if nrow(vis) > 0
        pts1 = Point2f[]; col1 = RGBAf[]
        pts2 = Point2f[]; col2 = RGBAf[]
        for r in eachrow(vis)
            h = relheight(r.I)
            (cr, cg, cb) = wl2rgb(r.λ)
            a = 0.15 + 0.85 * h
            push!(pts1, Point2f(r.λ, 0), Point2f(r.λ, 1))
            push!(col1, RGBAf(cr, cg, cb, a), RGBAf(cr, cg, cb, a))
            push!(pts2, Point2f(r.λ, 0), Point2f(r.λ, h))
            push!(col2, RGBAf(cr, cg, cb, 1), RGBAf(cr, cg, cb, 1))
        end
        linesegments!(ax1, pts1; color = col1, linewidth = 1.5)
        linesegments!(ax2, pts2; color = col2, linewidth = 1.5)
    end
    if show_exp && elem == 11
        for λe in exp_lines
            (λlo <= λe <= λhi) || continue
            vlines!(ax1, [λe]; color = (:white, 0.45), linestyle = :dash, linewidth = 1)
            vlines!(ax2, [λe]; color = (:white, 0.45), linestyle = :dash, linewidth = 1)
        end
    end
    xlims!(ax1, λlo, λhi); xlims!(ax2, λlo, λhi)
    ylims!(ax1, 0, 1);     ylims!(ax2, 0, 1.05)
    fig
end

# ╔═╡ 95b49730-1baf-4ce3-8c3f-836672af2ffb
md"""
### Najsilniejsze linie w zakresie
"""

# ╔═╡ 72fe7767-8ba6-48df-a1de-84193032cdc4
begin
    top = first(sort(vis, :I; rev = true), 15)
    select(top, :symbol, :ion,
           :λ => ByRow(x -> round(x; digits = 4)) => :λ_nm,
           :I => ByRow(x -> round(x; sigdigits = 3)) => :I_użyte, :n_linii, :intensity, :Aki, :conf_i, :term_i, :J_i, :conf_k, :term_k, :J_k)
end

# ╔═╡ b8693e3a-dadb-4b3d-a93c-9900ccf2ab45
md"""
### Porównanie ze sprawozdaniem (Na, spektroskop pryzmatyczny)
Dla każdej linii z tabeli szukam w NIST najsilniejszej linii Na w oknie ±2 nm (długości fal w powietrzu).
Kolumna **dekad** mówi, ile rzędów wielkości poniżej najsilniejszej linii leży ta linia, czyli jaki **próg widoczności** jest potrzebny, żeby ją zobaczyć.
"""

# ╔═╡ 814de1f5-21b8-41e5-bc79-a2bdb96ed9b0
begin
    # λ [nm] z tabeli ze sprawozdania (x to położenie na skali spektroskopu)
    exp_lines = [653, 616, 589, 568, 516, 507, 441, 433, 429, 423]

    function match_exp(sub, λe, lmax; win = 2.0)
        λa = vac2air.(sub.wl_vac_nm)
        idx = findall(x -> abs(x - λe) <= win, λa)
        isempty(idx) && return (λ_exp = λe, λ_NIST = missing, Δλ = missing,
                                intensity = missing, dekad = missing, linii_w_oknie = 0)
        j = idx[argmax(sub.I[idx])]
        return (λ_exp = λe, λ_NIST = round(λa[j]; digits = 3), Δλ = round(λa[j] - λe; digits = 2),
                intensity = round(sub.I[j]; sigdigits = 3), dekad = round(lmax - log10(sub.I[j]); digits = 2),
                linii_w_oknie = length(idx))
    end
end

# ╔═╡ 3b785249-038d-4d07-bc86-01fcccd05b87
if elem == 11
    DataFrame([match_exp(sub, λe, lmax) for λe in exp_lines])
else
    md"Tabela ze sprawozdania dotyczy sodu (Z = 11). Wybierz Na, żeby zobaczyć porównanie."
end

# ╔═╡ Cell order:
# ╟─81eb466a-8faa-420b-a51e-36781af7868b
# ╠═9735532e-8621-4f11-bec1-f557068e179e
# ╠═b3c367fa-571f-4406-b32b-e3bfc541f3f0
# ╠═9b1b3791-b66e-4828-9a43-f4d2b9aa7ae7
# ╟─7b9c7753-081d-4d70-ab3c-f20def9c4a5d
# ╠═79cead40-7703-4893-937f-a1b0b4e0b070
# ╟─fc11e62f-5c5a-40cc-ae64-71d9f7117dbd
# ╠═533e7dde-9e80-4606-abc8-5747a5d7d29c
# ╟─95b49730-1baf-4ce3-8c3f-836672af2ffb
# ╠═72fe7767-8ba6-48df-a1de-84193032cdc4
# ╟─b8693e3a-dadb-4b3d-a93c-9900ccf2ab45
# ╠═814de1f5-21b8-41e5-bc79-a2bdb96ed9b0
# ╠═3b785249-038d-4d07-bc86-01fcccd05b87