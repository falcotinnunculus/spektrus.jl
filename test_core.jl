# test_core.jl - szybki test rdzenia na prawdziwych danych.
#   julia --project=. .\test_core.jl
# Wymaga spectrum.csv (po prepare_lines.jl) i Spektrus.jl w tym samym folderze.

include("Spektrus.jl")
using .Spektrus, DataFrames

ok = Ref(true)
function check(name, cond)
    println(rpad(name, 70), cond ? "OK" : "BŁĄD")
    cond || (ok[] = false)
end
sumI(lines, λ0; win = 0.5) = sum(lines.I[abs.(lines.λ .- λ0) .<= win]; init = 0.0)

df = load_data()
println("Wczytano $(nrow(df)) linii\n")

check("vac2air: H-α próżnia 656.4696 → powietrze 656.288", abs(vac2air(656.4696) - 656.2883) < 0.001)
check("vac2air: poniżej 200 nm bez zmian", vac2air(150.0) == 150.0)

# --- sód ---
na = prepare_lines(df, 11; source = "lte", T = 6000.0)
check("Na: tryb LTE aktywny", na.use_lte)
vis = na.lines[(na.lines.λ .>= 380) .& (na.lines.λ .<= 780), :]
top = vis[argmax(vis.I), :]
check("Na: najsilniejsza linia widzialna to D2 (588.995 nm)", abs(top.λ - 588.995) < 0.02)
Iref = reference_intensity(na.lines, 380, 780)
check("Na: linia odniesienia (visible) = D2", Iref == top.I)
check("Na: wysokość linii odniesienia = 1", line_height(Iref, Iref, 4.0) == 1.0)
check("Na: linia 4 dekady niżej = minimalna wysokość 0.04", line_height(Iref * 1e-4, Iref, 4.0) == 0.04)

w1 = cluster_window(na.lines, 588.5, 590.0, 1.0)
w2 = cluster_window(na.lines, 588.5, 590.0, 0.2)
c1 = argmax(w1.I)
check("Na: dublet D zlewa się przy 1 nm (suma > pojedyncza linia)", maximum(w1.I) > maximum(w2.I))
check("Na: zlana kreska leży między D2 a D1 (588.995–589.592)", 588.995 < w1.λ[c1] < 589.592)
check("Na: przy 0.2 nm D1 i D2 są osobno", length(w2.λ) >= 2)

# --- wodór ---
h = prepare_lines(df, 1; source = "lte", T = 6000.0)
check("H: tryb LTE aktywny (są składowe z J i Aki)", h.use_lte)
Ha, Hb, Hg = sumI(h.lines, 656.28), sumI(h.lines, 486.13), sumI(h.lines, 434.05)
check("H: dekrement Balmera Hα > Hβ > Hγ > 0", Ha > Hb > Hg > 0)
check("H: Hα/Hβ w LTE 6000 K między 4 a 15 (wyszło $(round(Ha / Hb; digits = 1)))", 4 < Ha / Hb < 15)

# --- powrót do NIST i pierwiastki specjalne ---
for (Z, nm) in [(90, "Th"), (92, "U")]
    r = prepare_lines(df, Z; source = "lte")
    check("$nm: brak Aki → natężenia NIST (use_lte == false, są linie)", !r.use_lte && nrow(r.lines) > 0)
end
fr = prepare_lines(df, 87; source = "lte")
check("Fr: LTE działa mimo braku natężeń NIST", fr.use_lte && nrow(fr.lines) > 0)
check("Fr: tryb NIST pusty (brak natężeń)", nrow(prepare_lines(df, 87; source = "nist").lines) == 0)

# --- wszystkie pierwiastki dają coś do narysowania ---
els = element_list(df)
empty_els = [r.Z for r in eachrow(els) if nrow(prepare_lines(df, r.Z; source = "lte").lines) == 0]
check("Każdy z $(nrow(els)) pierwiastków ma linie do narysowania (LTE lub NIST)", isempty(empty_els))
isempty(empty_els) || println("   puste: ", empty_els)

# --- opis linii ---
println("\nPrzykład opisu: ", describe_line(na.lines, findfirst(==(top.λ), na.lines.λ)))

println(ok[] ? "\nWszystkie testy przeszły." : "\nSą błędy.")
exit(ok[] ? 0 : 1)