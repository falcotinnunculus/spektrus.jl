# test_core.jl - szybki test rdzenia (Spektrus.jl) na prawdziwych danych; interfejsu nie testuje.
#   julia --project=. .\test_core.jl
# Wymaga spectrum.csv (po prepare_lines.jl) i Spektrus.jl w tym samym folderze. levels.csv nie jest potrzebny:
# testy wyłączają go (poza sekcją o levels.csv, która używa własnego, sztucznego pliku).

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

# testy poniżej (poza sekcją o levels.csv) liczymy na poziomach z linii, niezależnie od tego, czy levels.csv istnieje
Spektrus.LEVEL_DB[] = false

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

# --- stopień bez Aki w trybie LTE ---
let r = prepare_lines(df, 57; ions = [1, 3], source = "lte")
    check("La I+III w LTE: stopień III pominięty (brak Aki) i zgłoszony w `skipped`",
          r.use_lte && haskey(r.skipped, 3) && !any(==(3), r.lines.ion) && any(==(1), r.lines.ion))
    check("Na I w LTE: nic nie jest zgłoszone jako całkowicie pominięte", isempty(prepare_lines(df, 11; ions = [1], source = "lte").skipped))
end

# --- równowaga Saha ---
# wartości odniesienia policzone niezależnie w Pythonie (SAHA_C z CODATA = 4.8294e15)
fH = saha_fractions([2.0, 1.0], [109678.77], 10000.0, 1e14)
check("Saha: wodór, 10000 K, n_e = 1e14: f_I = 0.2281 (wyszło $(round(fH[1]; digits = 4)))", abs(fH[1] - 0.22805) < 1e-3)
check("Saha: udziały sumują się do 1", abs(sum(fH) - 1) < 1e-12)
fNa = saha_fractions([2.1028, 1.0], [41449.45], 6000.0, 1e15)
check("Saha: sód, 6000 K, n_e = 1e15: f_I = 0.0191 (wyszło $(round(fNa[1]; digits = 4)))", abs(fNa[1] - 0.01905) < 5e-4)
check("Saha: większe n_e cofa jonizację", saha_fractions([2.0, 1.0], [109678.77], 10000.0, 1e16)[1] > fH[1])

ion = ionization_table()
if ion === nothing
    println("(brak ionization_energies.csv: testy Sahy na danych pominięte)")
else
    check("Energia jonizacji Na I = 41449.45 cm⁻¹", abs(ion[(11, 1)] - 41449.45) < 0.1)
    na_s = prepare_lines(df, 11; source = "lte", T = 6000.0, saha = true, ne = 1e15)
    check("Na: Saha aktywna (są udziały)", na_s.fractions !== nothing)
    na_s.fractions !== nothing &&
        check("Na: przy 6000 K i n_e = 1e15 neutralnych jest ok. 2% (wyszło $(round(100 * na_s.fractions[1]; digits = 2))%)",
              0.01 < na_s.fractions[1] < 0.04)
    rr = na_s.lines.I ./ na.lines.I
    check("Na: Saha skaluje wszystkie linie neutralne tym samym czynnikiem", maximum(rr) / minimum(rr) < 1 + 1e-9)
    hs = prepare_lines(df, 1; source = "lte", T = 10000.0, saha = true, ne = 1e14)
    check("H: Saha uwzględnia proton (f_I ≈ 0.23 przy 10000 K, n_e = 1e14)",
          hs.fractions !== nothing && abs(hs.fractions[1] - 0.228) < 0.01)
end

# --- jednostki i nazwy ---
check("jednostki: 500 nm = 20000 cm⁻¹", unit_value(500.0, "cm⁻¹") ≈ 20000.0)
check("jednostki: 500 nm = 2.4797 eV", abs(unit_value(500.0, "eV") - 2.479684) < 1e-5)
check("jednostki: 500 nm = 5000 Å = 0.5 µm", unit_value(500.0, "Å") ≈ 5000.0 && unit_value(500.0, "µm") ≈ 0.5)
check("jednostki: konwersja w obie strony",
      all(abs(unit_to_nm(unit_value(633.2, u), u) - 633.2) < 1e-9 for u in ("nm", "Å", "µm", "cm⁻¹", "eV")))
check("nazwy: Na = Sód, Fe = Żelazo, Es = Einstein",
      element_name(11) == "Sód" && element_name(26) == "Żelazo" && element_name(99) == "Einstein")

check("nazwy EN: Na = Sodium, Fe = Iron, Es = Einsteinium",
      element_name(11; lang = "en") == "Sodium" && element_name(26; lang = "en") == "Iron" && element_name(99; lang = "en") == "Einsteinium")
let i = findfirst(==(top.λ), na.lines.λ)
    check("opis linii: domyślnie po polsku, lang = en po angielsku (zlane linie)",
          occursin("zlane linie", describe_line(na.lines, i; n = 2)) && occursin("merged lines", describe_line(na.lines, i; n = 2, lang = "en")))
    check("opis linii: tryb energii, EN = vacuum, PL = próżnia",
          occursin("próżnia", describe_line(na.lines, i; unit = "eV")) && occursin("vacuum", describe_line(na.lines, i; unit = "eV", lang = "en")))
end

# --- pełne listy poziomów (levels.csv) ---
let tmp = tempname() * ".csv"
    write(tmp, "Z,ion,E_cm1,g\n11,1,0.0,2\n11,1,16956.17,2\n11,1,16973.37,4\n11,1,60000.0,2\n")   # ostatni leży powyżej chi
    db = load_levels(tmp)
    check("levels.csv: wczytanie (Na I, 4 poziomy)", db !== nothing && length(db[(11, 1)]) == 4)
    base = level_table(df, 11, 1, 41449.45)                    # tryb poziomów z linii (LEVEL_DB = false)
    Spektrus.LEVEL_DB[] = db
    lv = level_table(df, 11, 1, 41449.45)                      # cache rozróżnia bazy poziomów, więc to nowy wynik
    check("level_table: poziom podstawowy z NIST (g = 2), bez sztucznego (0.0, 1)", (0.0, 2) in lv && !((0.0, 1) in lv))
    check("level_table: poziomy powyżej energii jonizacji odcięte", all(p -> p[1] < 41449.45, lv))
    check("level_table: poziom z NIST nie jest liczony dwa razy (zgodny z poziomem z linii)",
          count(p -> abs(p[1] - 16973.37) < 0.01, lv) == 1)
    check("level_table: NIST + poziomy z linii (nie mniej niż same linie)", length(lv) >= length(base) - 1 && length(lv) > 3)
    Spektrus.LEVEL_DB[] = false
    check("bez levels.csv: powrót do poziomów z linii", level_table(df, 11, 1, 41449.45) == base)
    rm(tmp)
end

# --- tabela linii i adres strony ---
let r = prepare_lines(df, 11; ions = [1], source = "lte")
    Iref = reference_intensity(r.lines, 380.0, 780.0)
    n = size(r.lines, 1)
    tb = line_table(r.lines, 1:n, Iref, "nm")
    check("line_table: wiersz na każdą linię, od najsilniejszej", size(tb, 1) == n && issorted(tb.intensity_percent; rev = true))
    d2 = tb[abs.(tb.wavelength_nm .- 588.995) .< 0.01, :]
    check("line_table: D2 ma ok. 100% (linia odniesienia) i poziom dolny 3s",
          !isempty(d2) && maximum(d2.intensity_percent) ≈ 100 && occursin("3s", d2.lower_level[1]))
    te = line_table(r.lines, 1:n, Iref, "eV")
    check("line_table: w eV wartość = hc/λ(próżnia)", all(abs.(te.axis_value .- 1239.841984 ./ te.wavelength_vac_nm) .< 1e-9))
end
let q = parse_query("/?z=26&x=380,780&lang=en&ions=1%2C3")
    check("parse_query: klucze i wartości", q["z"] == "26" && q["x"] == "380,780" && q["lang"] == "en" && q["ions"] == "1,3")
    check("parse_query: brak zapytania daje pusty słownik", isempty(parse_query("/")) && isempty(parse_query("/?")))
    check("snap: wartość z adresu trafia w krok suwaka", snap(1000:100:20000, 6049) == 6000 && snap(0.5:0.1:6.0, 4.04) ≈ 4.0)
end

# --- język z przeglądarki ---
check("język: pl-PL na pierwszym miejscu daje polski", pick_language("pl-PL,pl;q=0.9,en-US;q=0.8,en;q=0.7") == "pl")
check("język: en-US przed pl daje angielski", pick_language("en-US,en;q=0.9,pl;q=0.8") == "en")
check("język: brak polskiego i angielskiego daje angielski (domyślny)", pick_language("de-DE,de;q=0.9,fr;q=0.5") == "en")
check("język: polski dalej na liście, ale jedyny obsługiwany, daje polski", pick_language("de,pl;q=0.5") == "pl")
check("język: pusty nagłówek, gwiazdka i q=0 dają domyślny", pick_language("") == "en" && pick_language("*") == "en" && pick_language("pl;q=0, en;q=0.5") == "en")

# --- opis linii ---
println("\nPrzykład opisu: ", describe_line(na.lines, findfirst(==(top.λ), na.lines.λ)))

println(ok[] ? "\nWszystkie testy przeszły." : "\nSą błędy.")
exit(ok[] ? 0 : 1)