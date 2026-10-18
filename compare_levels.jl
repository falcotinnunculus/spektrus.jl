# compare_levels.jl
# Porownuje funkcje podzialu U_z i udzialy stopni jonizacji N_z/N liczone dwoma sposobami:
#   (a) tylko z poziomow wystepujacych w pobranych liniach (stary tryb),
#   (b) z pelnych list poziomow NIST (levels.csv) uzupelnionych poziomami z linii.
#
# Wymaga spectrum.csv, ionization_energies.csv i levels.csv (po fetch_nist_levels.jl).
#   julia --project=. compare_levels.jl            # domyslnie Na (11) i Fe (26)
#   julia --project=. compare_levels.jl 6 26 57    # wybrane pierwiastki (Z)

include("Spektrus.jl")
using .Spektrus, DataFrames, Printf

df = load_data()
db = load_levels()
db === nothing && error("Brak levels.csv - najpierw uruchom fetch_nist_levels.jl")
zs  = isempty(ARGS) ? [11, 26] : parse.(Int, ARGS)
ne  = 1e15
ROM = ["I", "II", "III", "IV", "V"]

for Z in zs
    println("\n=== Z = $Z ($(element_name(Z))), n_e = $ne cm^-3 ===")
    for T in (6000.0, 10000.0, 15000.0)
        Spektrus.LEVEL_DB[] = false
        a = ionization_balance(df, Z, T, ne)
        Spektrus.LEVEL_DB[] = db                         # cache poziomow rozroznia tryby (klucz zawiera baze poziomow)
        b = ionization_balance(df, Z, T, ne)
        (a === nothing || b === nothing) && (println("T = $(Int(T)) K: brak danych"); continue)
        println("T = $(Int(T)) K   (liczba uwzglednionych stopni: linie $(length(a.U)), pelne $(length(b.U)))")
        println("  stopien |  U (linie)   U (pelne)   U pelne/linie |  udzial (linie)  udzial (pelne)")
        for z in 1:max(length(a.U), length(b.U))
            ua = z <= length(a.U) ? @sprintf("%10.3f", a.U[z]) : "         -"
            ub = z <= length(b.U) ? @sprintf("%10.3f", b.U[z]) : "         -"
            ra = z <= length(a.U) && z <= length(b.U) ? @sprintf("%12.2f", b.U[z] / a.U[z]) : "           -"
            fa = z <= length(a.f) ? @sprintf("%13.4f", a.f[z]) : "            -"
            fb = z <= length(b.f) ? @sprintf("%13.4f", b.f[z]) : "            -"
            println("  ", rpad(ROM[z], 7), " | ", ua, " ", ub, " ", ra, "   | ", fa, "  ", fb)
        end
    end
end
Spektrus.LEVEL_DB[] = nothing