# fetch_nist_lines.jl
# Pobiera linie emisyjne pierwiastkow 1-99 (stopnie jonizacji I i II) z NIST ASD
# i sklada je w jeden plik lines.csv.
#
# Instalacja pakietow (raz):
#   julia -e 'using Pkg; Pkg.add(["HTTP", "CSV", "DataFrames", "ProgressMeter"])'
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

using HTTP, CSV, DataFrames, ProgressMeter

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

function parse_lines(body, Z, sym)
    rows = DataFrame(
        Z = Int[], symbol = String[], ion = Int[],
        wl_nm = Float64[], medium = String[], kind = String[],
        intensity = Union{Missing,Float64}[],
        Aki = String[], conf_i = String[], term_i = String[], J_i = String[],
        conf_k = String[], term_k = String[], J_k = String[],
    )
    (isempty(body) || occursin("<html", lowercase(first(body, 300)))) && return rows

    lines = split(body, '\n')
    hi = findfirst(l -> occursin("obs_wl", l) || occursin("ritz_wl", l), lines)
    hi === nothing && return rows            # np. "No lines are available"

    header = [clean(h) for h in split(lines[hi], '\t')]
    findcol(prefix) = findfirst(h -> startswith(h, prefix), header)

    c_obs   = findcol("obs_wl")
    c_ritz  = findcol("ritz_wl")
    c_int   = findcol("intens")
    c_sp    = findcol("sp_num")
    extra   = [(name, key, findfirst(==(name), header)) for (name, key) in KEEP_TEXT]

    for l in lines[hi+1:end]
        isempty(strip(l)) && continue
        cells = [clean(c) for c in split(l, '\t'; keepempty = true)]
        length(cells) < length(header) && append!(cells, fill("", length(header) - length(cells)))
        startswith(cells[1], "element") && continue     # powtorzony naglowek

        wl = c_obs  === nothing ? missing : num(cells[c_obs])
        kind, cidx = "obs", c_obs
        if ismissing(wl)
            wl = c_ritz === nothing ? missing : num(cells[c_ritz])
            kind, cidx = "ritz", c_ritz
        end
        ismissing(wl) && continue

        medium = occursin("vac", header[cidx]) ? "vac" : "air"
        ion = c_sp === nothing ? 1 : something(tryparse(Int, cells[c_sp]), 1)
        inten = c_int === nothing ? missing : num(cells[c_int])

        get_extra(key) = begin
            for (_, k, ci) in extra
                k == key && ci !== nothing && return cells[ci]
            end
            ""
        end

        push!(rows, (Z, sym, ion, wl, medium, kind, inten,
                     get_extra(:Aki), get_extra(:conf_i), get_extra(:term_i), get_extra(:J_i),
                     get_extra(:conf_k), get_extra(:term_k), get_extra(:J_k)))
    end
    return rows
end

# --- main -------------------------------------------------------------------

function main()
    zmin = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
    zmax = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 99
    mkpath(RAWDIR)

    frames = DataFrame[]
    total_lines = 0
    n = zmax - zmin + 1
    prog = Progress(n;
        desc = "Pobieram z NIST: ",
        dt   = 0.5,
        showspeed = true,
        output = stderr)

    for Z in zmin:zmax
        sym = SYMBOLS[Z]
        cached = isfile(joinpath(RAWDIR, "$(lpad(Z, 2, '0'))_$(sym).txt"))
        next!(prog; showvalues = [
            (:Z, Z),
            (:symbol, sym),
            (:cached, cached),
            (:lines_so_far, total_lines),
        ])
        body = fetch_raw(Z, sym)
        df = parse_lines(body, Z, sym)
        total_lines += nrow(df)
        nrow(df) > 0 && push!(frames, df)
    end
    finish!(prog)

    if isempty(frames)
        @error "Nie udalo sie sparsowac zadnych linii. Zajrzyj do plikow w $RAWDIR/ - pokaz mi ich poczatek."
        return
    end
    all = reduce(vcat, frames)
    CSV.write(OUT, all)
    println("\nZapisano $(nrow(all)) linii do $OUT")
end

main()