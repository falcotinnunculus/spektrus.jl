# fetch_nist_lines.jl
# Pobiera linie emisyjne pierwiastkow 1-99 (stopnie jonizacji I i II) z NIST ASD
# i sklada je w jeden plik lines.csv.
#
# Instalacja pakietow (raz):
#   julia -e 'using Pkg; Pkg.add(["HTTP", "CSV", "DataFrames"])'
#
# Test na samym wodorze (zacznij od tego!):
#   julia fetch_nist_lines.jl 1 1
#   -> w lines.csv powinna byc m.in. linia H-alfa ok. 656.28 nm
#
# Pelne pobranie (ok. 99 zapytan, ~2-5 min, wynik jest cache'owany w raw/):
#   julia fetch_nist_lines.jl
#
# Wybrany zakres pierwiastkow:
#   julia fetch_nist_lines.jl 20 30
#
# Dane: NIST Atomic Spectra Database, https://physics.nist.gov/asd
# Cytowanie: Kramida, Ralchenko, Reader and NIST ASD Team, NIST Atomic Spectra Database.

using HTTP, CSV, DataFrames

const URL    = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"
const RAWDIR = "raw"
const OUT    = "lines.csv"
const LOW_NM = 100.0     # dolna granica zakresu [nm] (UV)
const UPP_NM = 2500.0    # gorna granica zakresu [nm] (bliska IR)
const PAUSE  = 1.0       # przerwa miedzy zapytaniami [s], badzmy grzeczni dla serwera

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

function build_query(sym)
    # Parametry formularza lines1.pl; format=3 -> tabulatory.
    # (Spisane z pamieci - jesli NIST zmienil formularz, tu trzeba poprawic.)
    return [
        "spectra"      => "$sym I-II",
        "low_w"        => string(LOW_NM),
        "upp_w"        => string(UPP_NM),
        "unit"         => "1",          # nm
        "submit"       => "Retrieve Data",
        "format"       => "3",          # tab-delimited
        "line_out"     => "0",          # wszystkie linie (filtrujemy sami)
        "en_unit"      => "0",
        "output"       => "0",          # calosc, bez stronicowania
        "bibrefs"      => "1",
        "page_size"    => "15",
        "show_obs_wl"  => "1",
        "show_calc_wl" => "1",
        "order_out"    => "0",
        "max_low_enrg" => "",
        "show_av"      => "2",          # proznia <200 nm, powietrze 200-2000 nm, proznia >2000 nm
        "max_upp_enrg" => "",
        "tsb_value"    => "0",
        "min_str"      => "",
        "A_out"        => "0",
        "f_out"        => "on",
        "intens_out"   => "on",
        "max_str"      => "",
        "allowed_out"  => "1",
        "forbid_out"   => "1",
        "min_accur"    => "",
        "min_intens"   => "",
        "conf_out"     => "on",
        "term_out"     => "on",
        "enrg_out"     => "on",
        "J_out"        => "on",
    ]
end

function fetch_raw(Z, sym)
    path = joinpath(RAWDIR, "$(lpad(Z, 2, '0'))_$(sym).txt")
    isfile(path) && return read(path, String)
    local body = ""
    for attempt in 1:3
        try
            resp = HTTP.get(URL; query = build_query(sym), readtimeout = 120, retry = false)
            body = String(resp.body)
            break
        catch e
            @warn "Blad zapytania dla $sym (proba $attempt)" exception = e
            sleep(3attempt)
        end
    end
    isempty(body) && return ""
    write(path, body)
    sleep(PAUSE)
    return body
end

# --- parsowanie -------------------------------------------------------------

const KEEP_TEXT = [
    "Aki(s^-1)" => :Aki,
    "conf_i" => :conf_i, "term_i" => :term_i, "J_i" => :J_i,
    "conf_k" => :conf_k, "term_k" => :term_k, "J_k" => :J_k,
]

function read_header(line)
    header = [clean(h) for h in split(line, '\t')]
    findcol(prefix) = findfirst(h -> startswith(h, prefix), header)
    return (header = header,
            obs = findcol("obs_wl"), ritz = findcol("ritz_wl"),
            int = findcol("intens"), sp = findcol("sp_num"),
            extra = Dict(key => findfirst(==(name), header) for (name, key) in KEEP_TEXT))
end

# NIST powtarza wiersz naglowka przy zmianie osrodka (proznia/powietrze),
# wiec trzeba go czytac od nowa za kazdym razem.
is_header(l) = occursin(r"^(element|obs_|ritz_)", l)

function parse_lines(body, Z, sym)
    rows = DataFrame(
        Z = Int[], symbol = String[], ion = Int[],
        wl_nm = Float64[], medium = String[], kind = String[],
        intensity = Union{Missing,Float64}[],
        Aki = String[], conf_i = String[], term_i = String[], J_i = String[],
        conf_k = String[], term_k = String[], J_k = String[],
    )
    sp_missing = false
    (isempty(body) || occursin("<html", lowercase(first(body, 300)))) && return rows, sp_missing

    H = nothing
    for l in split(body, '\n')
        isempty(strip(l)) && continue
        if is_header(l)
            H = read_header(l)
            if H.obs === nothing && H.ritz === nothing
                @warn "Naglowek bez kolumn z dlugoscia fali (inna jednostka?) - pomijam te wiersze" sym header = H.header
                H = nothing
            end
            continue
        end
        H === nothing && continue

        cells = [clean(c) for c in split(l, '\t'; keepempty = true)]
        n = length(H.header)
        length(cells) < n && append!(cells, fill("", n - length(cells)))

        wl = H.obs === nothing ? missing : num(cells[H.obs])
        kind, cidx = "obs", H.obs
        if ismissing(wl)
            wl = H.ritz === nothing ? missing : num(cells[H.ritz])
            kind, cidx = "ritz", H.ritz
        end
        ismissing(wl) && continue

        medium = occursin("vac", H.header[cidx]) ? "vac" : "air"
        if H.sp === nothing
            sp_missing = true               # tylko jeden gatunek w odpowiedzi
            ion = 1
        else
            ion = something(tryparse(Int, cells[H.sp]), 1)
        end
        inten = H.int === nothing ? missing : num(cells[H.int])
        extra(key) = (ci = H.extra[key]; ci === nothing ? "" : cells[ci])

        push!(rows, (Z, sym, ion, wl, medium, kind, inten,
                     extra(:Aki), extra(:conf_i), extra(:term_i), extra(:J_i),
                     extra(:conf_k), extra(:term_k), extra(:J_k)))
    end
    return rows, sp_missing
end

# --- main -------------------------------------------------------------------

function main()
    zmin = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
    zmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 99
    mkpath(RAWDIR)

    frames = DataFrame[]
    for Z in zmin:zmax
        sym = SYMBOLS[Z]
        body = fetch_raw(Z, sym)
        df, sp_missing = parse_lines(body, Z, sym)
        println(rpad("Z=$Z $sym", 12), nrow(df), " linii", sp_missing ? "  (brak sp_num -> ion=1)" : "")
        nrow(df) > 0 && push!(frames, df)
    end

    if isempty(frames)
        @error "Nie udalo sie sparsowac zadnych linii. Zajrzyj do plikow w $RAWDIR/ - pokaz mi ich poczatek."
        return
    end
    all = reduce(vcat, frames)
    CSV.write(OUT, all)
    println("\nZapisano $(nrow(all)) linii do $OUT")
    println(combine(groupby(all, [:medium, :kind]), nrow => :n,
                    :wl_nm => minimum => :min_nm, :wl_nm => maximum => :max_nm))
end

main()