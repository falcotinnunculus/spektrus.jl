# Spektrus.jl - rdzeń aplikacji (bez interfejsu): dane, fizyka, łączenie linii.
# Używany przez test_core.jl i (w następnym kroku) przez aplikację w przeglądarce.
#
#   include("Spektrus.jl"); using .Spektrus

module Spektrus

using CSV, DataFrames

export load_data, element_list, vac2air, wl2rgb, cluster_ids, prepare_lines,
       cluster_window, reference_intensity, visible_mask, line_height, describe_line,
       load_ionization, ionization_table, level_table, partition_function,
       saha_fractions, ionization_balance

const C2  = 1.438777            # hc/k [cm·K]
const VIS = (380.0, 780.0)      # zakres widzialny [nm]
const SAHA_C    = 4.8294e15     # 2 (2π m_e k / h²)^{3/2} [cm⁻³ K^{-3/2}] (CODATA)
const MAX_STAGE = 5             # najwyższy rozważany stopień jonizacji (V)

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

"Moment pędu J z tekstu: \"3/2\" → 1.5, \"1\" → 1.0; `missing`, gdy nie da się odczytać."
function parse_J(s)
    ismissing(s) && return missing
    t = strip(s)
    if occursin('/', t)
        a = split(t, '/')
        length(a) == 2 || return missing
        n = tryparse(Float64, a[1]); d = tryparse(Float64, a[2])
        (n === nothing || d === nothing || d == 0) && return missing
        return n / d
    end
    x = tryparse(Float64, t)
    return x === nothing ? missing : x
end

"Wczytuje spectrum.csv (wynik prepare_lines.jl); dolicza `gi` (waga statystyczna poziomu dolnego), jeśli go brak."
function load_data(path = joinpath(@__DIR__, "spectrum.csv"))
    str_cols = Dict(c => String for c in
        [:symbol, :kind, :Aki, :conf_i, :term_i, :J_i, :conf_k, :term_k, :J_k])
    df = CSV.read(path, DataFrame; types = str_cols)
    if !(:gi in propertynames(df))
        df.gi = [(j = parse_J(x); ismissing(j) ? missing : round(Int, 2j + 1)) for x in df.J_i]
    end
    return df
end

"Lista pierwiastków obecnych w danych: kolumny Z, symbol, n."
element_list(df) = sort(combine(groupby(df, [:Z, :symbol]), nrow => :n), :Z)

# ---------------- równowaga jonizacji (Saha) ----------------

const IONIZATION = Ref{Any}(nothing)     # (Z, stopień) => energia jonizacji [cm⁻¹]
const LEVELS = Dict{Tuple{UInt,Int,Int},Vector{Tuple{Float64,Int}}}()

"""
Wczytuje ionization_energies.csv (kolumny Z, stage, chi_cm1): energia jonizacji stopnia `stage`
(stage = 1 to I → II) w cm⁻¹. Zwraca `nothing`, gdy pliku nie ma.
"""
function load_ionization(path = joinpath(@__DIR__, "ionization_energies.csv"))
    isfile(path) || return nothing
    t = CSV.read(path, DataFrame)
    return Dict{Tuple{Int,Int},Float64}((r.Z, r.stage) => r.chi_cm1 for r in eachrow(t))
end

function ionization_table()
    IONIZATION[] === nothing && (IONIZATION[] = load_ionization())
    return IONIZATION[]
end

"""
Poziomy (E [cm⁻¹], g) stopnia `z` pierwiastka `Z`, zebrane z dolnych i górnych poziomów pobranych linii
(tylko poziomy poniżej energii jonizacji `chi_cm1`). Poziom podstawowy bywa poza zakresem pobranych linii
(np. rezonansowe poniżej 100 nm), więc gdy go brak, dokładamy (0, g = 1): przybliżenie dobre dla
atomów i jonów o zamkniętych powłokach, dla pozostałych może się mylić o mały czynnik.
"""
function level_table(df, Z, z, chi_cm1 = Inf)
    key = (objectid(df), Z, z)
    return get!(LEVELS, key) do
        sub = df[(df.Z .== Z) .& (df.ion .== z), :]
        seen = Set{Tuple{Float64,Int}}()
        for (E, g) in zip(sub.Ei_cm1, sub.gi)
            (ismissing(E) || ismissing(g)) || push!(seen, (Float64(E), Int(g)))
        end
        for (E, g) in zip(sub.Ek_cm1, sub.gk)
            (ismissing(E) || ismissing(g)) || push!(seen, (Float64(E), Int(g)))
        end
        if isempty(seen)
            Tuple{Float64,Int}[]
        else
            any(p -> p[1] == 0.0, seen) || push!(seen, (0.0, 1))
            [p for p in sort!(collect(seen)) if p[1] < chi_cm1]
        end
    end
end

"Funkcja podziału U(T) = Σ g·exp(−E/kT) z listy poziomów (E w cm⁻¹)."
partition_function(levels, T) = sum(g * exp(-C2 * E / T) for (E, g) in levels; init = 0.0)

"""
Udziały N_z/N_całk stopni jonizacji z równowagi Saha.
`U[z]` to funkcje podziału stopni 1..K, `chi[z]` to energia jonizacji z → z+1 [cm⁻¹] (K − 1 wartości),
`T` w K, `ne` w cm⁻³. N_{z+1}·n_e / N_z = SAHA_C·T^{3/2}·(U_{z+1}/U_z)·exp(−χ/kT).
"""
function saha_fractions(U, chi, T, ne)
    K = length(U)
    logw = zeros(K)
    for z in 1:K-1
        logR = log(SAHA_C) + 1.5 * log(T) + log(U[z+1] / U[z]) - log(ne) - C2 * chi[z] / T
        logw[z+1] = logw[z] + logR
    end
    w = exp.(logw .- maximum(logw))
    return w ./ sum(w)
end

"""
Równowaga jonizacji pierwiastka `Z` w temperaturze `T` i gęstości elektronowej `ne`: zwraca (f, U),
czyli udziały i funkcje podziału stopni I..K, gdzie K to najwyższy stopień, dla którego są dane
(maks. V) i znana energia jonizacji poprzedniego stopnia. Stopnie powyżej K są pomijane.
Zwraca `nothing`, gdy brak ionization_energies.csv albo kolumny `gi`.
"""
function ionization_balance(df, Z, T, ne)
    chi = ionization_table()
    (chi === nothing || !(:gi in propertynames(df))) && return nothing
    U = Float64[]; χ = Float64[]
    for z in 1:MAX_STAGE
        z > 1 && !haskey(chi, (Z, z - 1)) && break
        u = nothing
        if z == Z + 1
            u = 1.0                                   # goły rdzeń (np. H II = proton)
        else
            lv = level_table(df, Z, z, get(chi, (Z, z), Inf))
            isempty(lv) || (u = partition_function(lv, T))
        end
        u === nothing && break
        push!(U, u)
        z > 1 && push!(χ, chi[(Z, z - 1)])
    end
    isempty(U) && return nothing
    return (; f = saha_fractions(U, χ, T, ne), U = U)
end

"""
Linie jednego pierwiastka z natężeniami `I`, posortowane po długości fali `λ`.

source = "lte"  : I = g_k·A_ki·ν·exp(−E_k/kT) (emisja optycznie cienka w LTE);
                 gdy pierwiastek nie ma danych do modelu, wraca do natężeń NIST (use_lte = false)
source = "nist" : natężenia NIST (jednostki umowne), tylko linie, które je mają

saha = true (tylko w trybie LTE): natężenie uwzględnia udział stopnia jonizacji z równowagi Saha
oraz funkcję podziału, I = (N_z/N)·g_k·A_ki·ν·exp(−E_k/kT)/U_z, dla gęstości elektronowej `ne` [cm⁻³].
Dzięki temu jasności linii różnych stopni jonizacji są ze sobą porównywalne. Udziały wracają w
`fractions` (wektor dla stopni I..K albo `nothing`, gdy Saha się nie zastosowała).

Stopnie jonizacji: `ions` (np. [1], [2], [1, 2]) wybiera dokładnie te stopnie; bez `ions`
działa dawne `ion2` (false: tylko I, true: I i II). Pusty wybór daje pusty wynik.
"""
function prepare_lines(df, Z; ion2 = false, ions = nothing, medium = "air", source = "lte",
                       T = 6000.0, saha = false, ne = 1e15)
    stages = ions === nothing ? (ion2 ? (1, 2) : (1,)) : ions
    sub = df[(df.Z .== Z) .& in.(df.ion, Ref(stages)), :]
    sub.λ = medium == "air" ? vac2air.(sub.wl_vac_nm) : copy(sub.wl_vac_nm)
    use_lte = false
    fractions = nothing
    if source == "lte" && :Aki_s in propertynames(sub)
        lte = dropmissing(sub, [:Aki_s, :Ek_cm1, :gk])
        if nrow(lte) > 0
            use_lte = true
            sub = lte
            bal = saha ? ionization_balance(df, Z, T, ne) : nothing
            if bal === nothing
                sub.I = [g * A * (1e7 / w) * exp(-C2 * E / T)
                         for (g, A, w, E) in zip(sub.gk, sub.Aki_s, sub.wl_vac_nm, sub.Ek_cm1)]
            else
                fractions = bal.f
                K = length(bal.f)
                sub.I = [(z <= K ? bal.f[z] / bal.U[z] : 0.0) * g * A * (1e7 / w) * exp(-C2 * E / T)
                         for (z, g, A, w, E) in zip(sub.ion, sub.gk, sub.Aki_s, sub.wl_vac_nm, sub.Ek_cm1)]
            end
        end
    end
    if !use_lte
        sub = dropmissing(sub, :intensity)
        sub.I = Float64.(sub.intensity)
    end
    sub = sub[sub.I .> 0, :]
    sub = sub[sortperm(sub.λ), :]
    return (; lines = sub, use_lte, fractions)
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
    ion = ("I", "II", "III", "IV", "V")[clamp(r.ion, 1, 5)]
    lower = join(fmt.([r.conf_i, r.term_i, r.J_i]), " ")
    upper = join(fmt.([r.conf_k, r.term_k, r.J_k]), " ")
    extra = n > 1 ? "  (zlane linie: $n)" : ""
    return "$(r.symbol) $ion  $(round(r.λ; digits = 3)) nm   $lower → $upper   Aki = $(fmt(r.Aki)) s⁻¹$extra"
end

end # module