# Spektrus.jl - rdzeń aplikacji (bez interfejsu): dane, fizyka, łączenie linii.
# Używany przez test_core.jl i (w następnym kroku) przez aplikację w przeglądarce.
#
#   include("Spektrus.jl"); using .Spektrus

module Spektrus

using CSV, DataFrames

export load_data, element_list, vac2air, wl2rgb, cluster_ids, prepare_lines,
       cluster_window, reference_intensity, visible_mask, line_height, describe_line

const C2  = 1.438777            # hc/k [cm·K]
const VIS = (380.0, 780.0)      # zakres widzialny [nm]

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
        return (1.00, 0.50, 0.45)
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

"Numery klastrów dla posortowanych długości fal: linia należy do klastra, dopóki λ − λ_start <= res."
function cluster_ids(λ, res)
    ids = zeros(Int, length(λ))
    isempty(λ) && return ids
    b = 1
    start = λ[1]
    for i in eachindex(λ)
        if λ[i] - start > res
            b += 1
            start = λ[i]
        end
        ids[i] = b
    end
    return ids
end

"Wczytuje spectrum.csv (wynik prepare_lines.jl)."
function load_data(path = joinpath(@__DIR__, "spectrum.csv"))
    str_cols = Dict(c => String for c in
        [:symbol, :kind, :Aki, :conf_i, :term_i, :J_i, :conf_k, :term_k, :J_k])
    return CSV.read(path, DataFrame; types = str_cols)
end

"Lista pierwiastków obecnych w danych: kolumny Z, symbol, n."
element_list(df) = sort(combine(groupby(df, [:Z, :symbol]), nrow => :n), :Z)

"""
Linie jednego pierwiastka z natężeniami `I`, posortowane po długości fali `λ`.

source = "lte"  : I = g_k·A_ki·ν·exp(−E_k/kT) (emisja optycznie cienka w LTE);
                 gdy pierwiastek nie ma danych do modelu, wraca do natężeń NIST (use_lte = false)
source = "nist" : natężenia NIST (jednostki umowne), tylko linie, które je mają

Stopnie jonizacji: `ions` (np. [1], [2], [1, 2]) wybiera dokładnie te stopnie; bez `ions`
działa dawne `ion2` (false: tylko I, true: I i II). Pusty wybór daje pusty wynik.
"""
function prepare_lines(df, Z; ion2 = false, ions = nothing, medium = "air", source = "lte", T = 6000.0)
    stages = ions === nothing ? (ion2 ? (1, 2) : (1,)) : ions
    sub = df[(df.Z .== Z) .& in.(df.ion, Ref(stages)), :]
    sub.λ = medium == "air" ? vac2air.(sub.wl_vac_nm) : copy(sub.wl_vac_nm)
    use_lte = false
    if source == "lte" && :Aki_s in propertynames(sub)
        lte = dropmissing(sub, [:Aki_s, :Ek_cm1, :gk])
        if nrow(lte) > 0
            use_lte = true
            sub = lte
            sub.I = [g * A * (1e7 / w) * exp(-C2 * E / T)
                     for (g, A, w, E) in zip(sub.gk, sub.Aki_s, sub.wl_vac_nm, sub.Ek_cm1)]
        end
    end
    if !use_lte
        sub = dropmissing(sub, :intensity)
        sub.I = Float64.(sub.intensity)
    end
    sub = sub[sub.I .> 0, :]
    sub = sub[sortperm(sub.λ), :]
    return (; lines = sub, use_lte)
end

"""
Zlewa linie z okna [λlo, λhi] leżące bliżej niż `res` [nm] (jak pryzmat lub piksel).
Natężenia są sumowane, położenie to średnia ważona natężeniem.
Zwraca (λ, I, n, rows): `rows` to indeksy najsilniejszej linii każdego klastra w `lines`.
"""
function cluster_window(lines, λlo, λhi, res)
    idx = findall(x -> λlo <= x <= λhi, lines.λ)
    λ = lines.λ[idx]
    I = lines.I[idx]
    ids = cluster_ids(λ, res)
    nc = isempty(ids) ? 0 : maximum(ids)
    Ic = zeros(nc); λw = zeros(nc); n = zeros(Int, nc); best = zeros(Int, nc)
    for j in eachindex(ids)
        c = ids[j]
        Ic[c] += I[j]
        λw[c] += λ[j] * I[j]
        n[c] += 1
        if best[c] == 0 || I[j] > I[best[c]]
            best[c] = j
        end
    end
    λw ./= Ic
    return (; λ = λw, I = Ic, n = n, rows = idx[best])
end

"""
Natężenie linii odniesienia dla skali jasności.
ref = "visible" (najsilniejsza linia 380–780 nm), "window" (w oknie λlo..λhi), "all" (cały pierwiastek).
"""
function reference_intensity(lines, λlo, λhi; ref = "visible")
    isempty(lines.I) && return 1.0
    sel = if ref == "all"
        trues(nrow(lines))
    elseif ref == "window"
        (lines.λ .>= λlo) .& (lines.λ .<= λhi)
    else
        (lines.λ .>= VIS[1]) .& (lines.λ .<= VIS[2])
    end
    any(sel) || (sel = trues(nrow(lines)))
    return maximum(lines.I[sel])
end

"Maska linii nie słabszych niż `decades` rzędów wielkości poniżej linii odniesienia."
visible_mask(I, Iref, decades) = log10.(I) .>= log10(Iref) - decades

"Wysokość kreski w [0.04, 1]: skala logarytmiczna (z progiem) albo liniowa."
function line_height(I, Iref, decades; scale = "log")
    scale == "log" ?
        clamp((log10(I) - (log10(Iref) - decades)) / decades, 0.04, 1.0) :
        min(I / Iref, 1.0)
end

fmt(x) = ismissing(x) ? "?" : string(x)

"Opis linii do podpowiedzi: pierwiastek, λ, przejście, Aki. `n` > 1 oznacza zlane linie."
function describe_line(lines, row; n = 1)
    r = lines[row, :]
    ion = ("I", "II", "III", "IV")[clamp(r.ion, 1, 4)]
    lower = join(fmt.([r.conf_i, r.term_i, r.J_i]), " ")
    upper = join(fmt.([r.conf_k, r.term_k, r.J_k]), " ")
    extra = n > 1 ? "  (zlane linie: $n)" : ""
    return "$(r.symbol) $ion  $(round(r.λ; digits = 3)) nm   $lower → $upper   Aki = $(fmt(r.Aki)) s⁻¹$extra"
end

end # module