# prepare_lines.jl
# Przygotowuje lines.csv (z fetch_nist_lines.jl) do spectrum.csv dla notebooka/aplikacji.
#
#   julia --project=. .\prepare_lines.jl
#
# Co robi:
#   1. sprowadza wszystkie długości fal do próżni (NIST daje powietrze w zakresie 200-2000 nm;
#      odwrotna konwersja Edlena/Mortona), przełącznik próżnia/powietrze jest w aplikacji
#   2. liczy kolumny liczbowe do modelu LTE: Aki_s [1/s], Ei_cm1, Ek_cm1 [cm^-1], gk = 2J+1
#   3. NIE wyrzuca linii bez natężenia NIST - model LTE ich nie potrzebuje, a wiele linii
#      (np. składowe struktury subtelnej wodoru) ma Aki i poziomy, ale nie ma natężenia.
#      Tryb "NIST" w aplikacji filtruje linie bez natężenia sam.
#   4. nie scala duplikatów długości fal - to robi aplikacja (rozdzielczość)

using CSV, DataFrames

const IN  = "lines.csv"
const OUT = "spectrum.csv"

# współczynnik załamania standardowego powietrza dla długości fali w próżni [nm]
# (bez zabezpieczenia zakresu: wywołujemy to tylko dla linii, które NIST podał w powietrzu)
function n_air(λvac)
    s2 = (1000 / λvac)^2
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
end

function air2vac(λair)
    λv = λair
    for _ in 1:5
        λv = λair * n_air(λv)               # iteracja, zbiega w 2-3 krokach
    end
    return λv
end

function num(s)
    ismissing(s) && return missing
    m = match(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", s)
    v = m === nothing ? nothing : tryparse(Float64, m.match)
    return v === nothing ? missing : v
end

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

str_cols = Dict(c => String for c in
    [:Aki, :conf_i, :term_i, :J_i, :conf_k, :term_k, :J_k, :Ei, :Ek, :medium, :kind, :symbol])
df = CSV.read(IN, DataFrame; types = str_cols)

println("Rozkład medium w lines.csv:")
println(combine(groupby(df, :medium), nrow => :n))

df.wl_nm = [m == "air" ? air2vac(w) : w for (w, m) in zip(df.wl_nm, df.medium)]
rename!(df, :wl_nm => :wl_vac_nm)

df.intensity = [(ismissing(x) || x <= 0) ? missing : x for x in df.intensity]
df.log10_int = [ismissing(x) ? missing : log10(x) for x in df.intensity]
df.Aki_s     = [num(x) for x in df.Aki]
df.Ei_cm1    = [num(x) for x in df.Ei]
df.Ek_cm1    = [num(x) for x in df.Ek]
df.gk        = [(j = parse_J(x); ismissing(j) ? missing : round(Int, 2j + 1)) for x in df.J_k]

# podsumowanie
df.has_int = .!ismissing.(df.intensity)
df.has_lte = completecases(df, [:Aki_s, :Ek_cm1, :gk])
cnt = combine(groupby(df, :Z), nrow => :wszystkie,
              :has_int => sum => :z_natezeniem_NIST, :has_lte => sum => :do_modelu_LTE)

select!(df, :Z, :symbol, :ion, :wl_vac_nm, :kind, :intensity, :log10_int,
        :Aki_s, :Ek_cm1, :Ei_cm1, :gk,
        :Aki, :conf_i, :term_i, :J_i, :conf_k, :term_k, :J_k)
sort!(df, [:Z, :wl_vac_nm])
CSV.write(OUT, df)

println("Zapisano $(nrow(df)) linii do $OUT")
println("Pierwiastków z liniami: $(nrow(cnt)) / 99")
no_lte = cnt.Z[cnt.do_modelu_LTE .== 0]
no_int = cnt.Z[cnt.z_natezeniem_NIST .== 0]
isempty(no_lte) || println("Bez linii nadających się do modelu LTE: ", join(no_lte, ", "))
isempty(no_int) || println("Bez natężeń NIST (tryb NIST pusty): ", join(no_int, ", "))
println("Najwięcej linii: ", first(sort(cnt, :wszystkie, rev = true), 5))