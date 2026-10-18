# fetch_nist_levels.jl
# Pobiera PELNE listy poziomow energetycznych pierwiastkow 1-99 (stopnie jonizacji I-V) z NIST ASD
# i sklada je w plik levels.csv (kolumny: Z, ion, E_cm1, g). Plik jest opcjonalny: gdy go nie ma,
# Spektrus.jl liczy funkcje podzialu z poziomow wystepujacych w pobranych liniach (przyblizenie).
#
# Instalacja pakietow: te same co dla fetch_nist_lines.jl (HTTP, CSV, DataFrames).
#
# Test na jednym pierwiastku (zacznij od tego!):
#   julia --project=. fetch_nist_levels.jl 26 26
#   -> w raw_levels/ pojawia sie 26_Fe_I.txt ... 26_Fe_V.txt, w levels.csv poziomy zelaza.
#   Jesli cos sie nie zgadza, pokaz poczatek np. raw_levels/26_Fe_I.txt (pierwsze ~20 linii).
#
# Pelne pobranie (ok. 480 zapytan, ~10 min, wynik jest cache'owany w raw_levels/):
#   julia --project=. fetch_nist_levels.jl
#
# Wybrany zakres pierwiastkow:
#   julia --project=. fetch_nist_levels.jl 20 30
#
# Dane: NIST Atomic Spectra Database, https://physics.nist.gov/asd
# Cytowanie: Kramida, Ralchenko, Reader and NIST ASD Team, NIST Atomic Spectra Database.

using HTTP, CSV, DataFrames

const URL     = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl"
const RAWDIR  = "raw_levels"
const OUT     = "levels.csv"
const MAX_ION = 5          # stopnie jonizacji I..V (jak w fetch_nist_lines.jl)
const PAUSE   = 1.0        # przerwa miedzy zapytaniami [s]
const ROMAN   = ["I", "II", "III", "IV", "V"]

const SYMBOLS = split("""
H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn
Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La
Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po
At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es
""")
@assert length(SYMBOLS) == 99

# --- pomocnicze -------------------------------------------------------------

clean(s) = strip(replace(String(s), r"[\"=]" => ""))

function num(s)
    m = match(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", s)
    m === nothing ? missing : parse(Float64, m.match)
end

# J jako tekst ("3/2", "2") -> waga statystyczna g = 2J + 1 (lub missing)
function g_from_J(s)
    t = strip(s)
    isempty(t) && return missing
    if occursin('/', t)
        a = split(t, '/')
        length(a) == 2 || return missing
        n = tryparse(Float64, strip(a[1])); d = tryparse(Float64, strip(a[2]))
        (n === nothing || d === nothing || d == 0) && return missing
        return round(Int, 2 * n / d + 1)
    end
    x = tryparse(Float64, t)
    return x === nothing ? missing : round(Int, 2x + 1)
end

function build_query(sym, ion)
    # Parametry formularza energy1.pl; format=3 -> tabulatory.
    # (Spisane z pamieci - jesli NIST zmienil formularz, tu trzeba poprawic.)
    return [
        "spectrum"          => "$sym $(ROMAN[ion])",
        "units"             => "0",          # cm-1
        "format"            => "3",          # tab-delimited
        "output"            => "0",          # calosc, bez stronicowania
        "page_size"         => "15",
        "multiplet_ordered" => "0",
        "conf_out"          => "on",
        "term_out"          => "on",
        "level_out"         => "on",
        "unc_out"           => "1",
        "j_out"             => "on",
        "g_out"             => "on",         # waga statystyczna
        "lande_out"         => "on",
        "perc_out"          => "on",
        "biblio"            => "on",
        "temp"              => "",
        "submit"            => "Retrieve Data",
    ]
end

function fetch_raw(Z, sym, ion)
    path = joinpath(RAWDIR, "$(lpad(Z, 2, '0'))_$(sym)_$(ROMAN[ion]).txt")
    isfile(path) && return read(path, String)
    local body = ""
    for attempt in 1:3
        try
            resp = HTTP.get(URL; query = build_query(sym, ion), readtimeout = 300, retry = false)
            body = String(resp.body)
            break
        catch e
            @warn "Blad zapytania dla $sym $(ROMAN[ion]) (proba $attempt)" exception = e
            sleep(3attempt)
        end
    end
    isempty(body) && return ""
    write(path, body)
    sleep(PAUSE)
    return body
end

# --- parsowanie -------------------------------------------------------------

# Zwraca (rows, limit_cm1, header). `limit_cm1` to energia jonizacji z wiersza "Limit" (lub missing).
function parse_levels(body, Z, ion)
    rows = DataFrame(Z = Int[], ion = Int[], E_cm1 = Float64[], g = Int[])
    limit = missing
    header = String[]
    (isempty(body) || occursin("<html", lowercase(first(body, 300)))) && return rows, limit, header

    cE = cg = cJ = cT = nothing
    for l in split(body, '\n')
        isempty(strip(l)) && continue
        cells = [clean(c) for c in split(l, '\t'; keepempty = true)]
        low = lowercase.(cells)
        # wiersz naglowka: ma kolumne "Level (cm-1)" (nazwy ustalamy po fragmentach, nie dokladnie)
        if any(c -> startswith(c, "level"), low)
            header = cells
            cE = findfirst(c -> startswith(c, "level") && occursin("cm", c), low)
            cE === nothing && (cE = findfirst(c -> startswith(c, "level"), low))
            cg = findfirst(==("g"), low)
            cJ = findfirst(==("j"), low)
            cT = findfirst(==("term"), low)
            continue
        end
        cE === nothing && continue
        length(cells) < cE && continue

        term = cT === nothing || length(cells) < cT ? "" : cells[cT]
        e = num(cells[cE])
        if lowercase(term) == "limit"
            ismissing(e) || (limit = ismissing(limit) ? e : min(limit, e))   # NIST podaje kilka limitów serii; energia jonizacji to najnizszy
            continue
        end
        ismissing(e) && continue                       # poziomy bez energii ("x+...")
        occursin(r"[xX]", cells[cE]) && continue       # energia wzgledna nieznana

        g = cg !== nothing && length(cells) >= cg ? num(cells[cg]) : missing
        g = ismissing(g) ? missing : round(Int, g)
        if ismissing(g) && cJ !== nothing && length(cells) >= cJ
            g = g_from_J(cells[cJ])
        end
        ismissing(g) && continue
        push!(rows, (Z, ion, e, g))
    end
    return rows, limit, header
end

# --- main -------------------------------------------------------------------

function main()
    zmin = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
    zmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 99
    mkpath(RAWDIR)

    chi = isfile("ionization_energies.csv") ?
          Dict((r.Z, r.stage) => r.chi_cm1 for r in eachrow(CSV.read("ionization_energies.csv", DataFrame))) :
          Dict{Tuple{Int,Int},Float64}()

    frames = DataFrame[]
    for Z in zmin:zmax
        sym = SYMBOLS[Z]
        for ion in 1:min(MAX_ION, Z)                  # stopnie powyzej Z nie maja poziomow (goly rdzen)
            body = fetch_raw(Z, sym, ion)
            df, limit, header = parse_levels(body, Z, ion)
            note = ""
            if nrow(df) == 0
                note = "  <- BRAK POZIOMOW" * (isempty(header) ? " (nie znaleziono naglowka z kolumna Level)" : " (naglowek: $(join(header, " | ")))")
            else
                ground = any(==(0.0), df.E_cm1) ? "" : "  (uwaga: brak poziomu 0.0)"
                ref = get(chi, (Z, ion), missing)
                cmp = (ismissing(limit) || ismissing(ref)) ? "" :
                      "  limit NIST $(round(limit; digits = 1)) vs tabela $(round(ref; digits = 1))" *
                      (abs(limit - ref) / ref > 0.01 ? "  <- ROZNICA > 1%" : "")
                few = nrow(df) <= 3 ? "  (malo poziomow w NIST: Spektrus.jl uzupelni je poziomami z linii)" : ""
                note = ground * cmp * few
                push!(frames, df)
            end
            println(rpad("Z=$Z $sym $(ROMAN[ion])", 14), lpad(nrow(df), 6), " poziomow", note)
        end
    end

    if isempty(frames)
        @error "Nie udalo sie sparsowac zadnych poziomow. Zajrzyj do plikow w $RAWDIR/ - pokaz mi ich poczatek."
        return
    end
    all = reduce(vcat, frames)
    CSV.write(OUT, all)
    println("\nZapisano $(nrow(all)) poziomow do $OUT")
end

main()