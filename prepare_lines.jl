# prepare_lines.jl
# Czyści lines.csv z fetch_nist_lines.jl do pliku spectrum.csv gotowego do rysowania.
#
#   julia prepare_lines.jl
#
# Co robi:
#   1. zostawia tylko linie obserwowane (kind == "obs") z dodatnim natężeniem
#   2. scala duplikaty (ta sama długość fali, ten sam jon) - bierze wiersz z największym natężeniem
#      (NIST podaje natężenie na wierszu multipletu, a składowe struktury subtelnej go nie mają)
#   3. dodaje log10_int, bo natężenia rozciągają się na wiele rzędów wielkości
#
# NIST podaje fale w próżni (<200 nm, >2000 nm) i w powietrzu (200-2000 nm).
# Tu sprowadzamy wszystko do próżni (odwrotna konwersja Edlena/Mortona dla linii "air"),
# a przełącznik próżnia/powietrze jest w aplikacji.

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

str_cols = Dict(c => String for c in
    [:Aki, :conf_i, :term_i, :J_i, :conf_k, :term_k, :J_k, :medium, :kind, :symbol])
df = CSV.read(IN, DataFrame; types = str_cols)

println("Rozkład medium w lines.csv:")
println(combine(groupby(df, :medium), nrow => :n))
df.wl_nm = [m == "air" ? air2vac(w) : w for (w, m) in zip(df.wl_nm, df.medium)]

filter!(r -> r.kind == "obs", df)
dropmissing!(df, :intensity)
filter!(r -> r.intensity > 0, df)
rename!(df, :wl_nm => :wl_vac_nm)

out = combine(groupby(df, [:Z, :ion, :wl_vac_nm])) do sdf
    i = argmax(sdf.intensity)
    (symbol = sdf.symbol[i], intensity = sdf.intensity[i], Aki = sdf.Aki[i],
     conf_i = sdf.conf_i[i], term_i = sdf.term_i[i], J_i = sdf.J_i[i],
     conf_k = sdf.conf_k[i], term_k = sdf.term_k[i], J_k = sdf.J_k[i])
end

out.log10_int = log10.(out.intensity)
select!(out, :Z, :symbol, :ion, :wl_vac_nm, :intensity, :log10_int,
        Not([:Z, :symbol, :ion, :wl_vac_nm, :intensity, :log10_int]))
sort!(out, [:Z, :wl_vac_nm])

CSV.write(OUT, out)

println("Zapisano $(nrow(out)) linii do $OUT")
counts = combine(groupby(out, :Z), nrow => :n)
println("Pierwiastków z liniami: $(nrow(counts)) / 99")
missingZ = setdiff(1:99, counts.Z)
isempty(missingZ) || println("Bez żadnych linii z natężeniem: ", join(missingZ, ", "))
println("Najwięcej linii: ", first(sort(counts, :n, rev = true), 5))