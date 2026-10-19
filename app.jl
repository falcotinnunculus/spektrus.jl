# app.jl - Spektruś.jl w przeglądarce (Bonito + WGLMakie).
#
# Instalacja (raz; w tym samym oknie PowerShella ogranicz równoległość, żeby nie zabrakło RAM-u):
#   $env:JULIA_NUM_PRECOMPILE_TASKS = 1
#   julia --project=.
#   ]add Bonito WGLMakie
#
# Uruchomienie:
#   julia --project=. .\app.jl
# i otwórz http://127.0.0.1:8080 (Ctrl+C kończy).
#
# Obsługa wykresu: kółko myszy = zoom, przeciągnięcie lewym = zaznaczenie obszaru do powiększenia,
# przeciągnięcie prawym = przesuwanie, przyciski UV/widzialne/IR/wszystko = zakres.

using Bonito, WGLMakie
include(joinpath(@__DIR__, "Spektrus.jl"))
using .Spektrus
using Logging
import CSV
import Base64

# Zwykłe kliknięcie na wykresie (bez przeciągania) daje zaznaczenie o zerowej szerokości; Makie ignoruje je
# i wypisuje ostrzeżenie "Rectangle zoom ignored". Jest nieszkodliwe, więc filtrujemy tylko tę jedną wiadomość.
struct DropZoomWarn{L<:AbstractLogger} <: AbstractLogger
    inner::L
end
Logging.min_enabled_level(l::DropZoomWarn) = Logging.min_enabled_level(l.inner)
Logging.shouldlog(l::DropZoomWarn, args...) = Logging.shouldlog(l.inner, args...)
Logging.catch_exceptions(l::DropZoomWarn) = Logging.catch_exceptions(l.inner)
function Logging.handle_message(l::DropZoomWarn, level, message, _module, group, id, file, line; kwargs...)
    occursin("Rectangle zoom ignored", string(message)) && return nothing
    return Logging.handle_message(l.inner, level, message, _module, group, id, file, line; kwargs...)
end
global_logger(DropZoomWarn(global_logger()))

# wykres wypełnia dostępne miejsce (rozmiar bierze z rodzica); starsze WGLMakie bez resize_to -> stały rozmiar
try
    WGLMakie.activate!(resize_to = :parent)
catch
    WGLMakie.activate!()
end

# ---- granice osi λ (zmień tutaj) ----
const X_MIN       = 0.0       # dolna granica [nm], skala liniowa
const X_MIN_LOG   = 10.0      # dolna granica [nm], skala logarytmiczna (log(0) nie istnieje)
const X_MAX       = 3000.0    # górna granica [nm]; dane NIST sięgają 2500 nm
const X_MIN_WIDTH = 0.002     # najwęższe okno [nm] w skali liniowej (w logarytmicznej: 0,002 dekady)
const TABLE_ROWS   = 200      # ile najsilniejszych linii pokazuje tabela (CSV zawiera wszystkie linie w widoku)
const LOAD_TIMING = false     # true: konsola pokazuje czasy kolejnych etapów budowania strony przy każdym wczytaniu
const RES_FLOOR   = 1e-4      # najmniejsza rozdzielczość zlewania linii [nm] (linie bliższe niż to uznajemy za tę samą)

# autor tej wersji (wyświetlany w oknie "O programie" obok Claude'a) - wpisz swoje imię i nazwisko
const AUTHOR = "Mateusz Winiarski"
const REPO_URL = "https://github.com/falcotinnunculus/spektrus.jl"

const DF       = load_data()
const ELEMENTS = element_list(DF)
const NAME_OPTS_PL = [element_name(Z; lang = "pl") for Z in ELEMENTS.Z]
const NAME_OPTS_EN = [element_name(Z; lang = "en") for Z in ELEMENTS.Z]
const SYM_OPTS  = collect(String, ELEMENTS.symbol)
# zakresy: (napis PL, napis EN, (od, do) [nm])
const RANGES   = [("UV", "UV", (100.0, 380.0)), ("widzialne", "visible", (380.0, 780.0)),
                  ("IR", "IR", (780.0, 2500.0)), ("wszystko", "all", (100.0, 2500.0))]

# kompilacja JIT dzieje się raz na proces Julii, więc przycisk "Rozgrzej" znika dla wszystkich nowych stron po pierwszym użyciu
const WARMED = Ref(false)

# język interfejsu: adres (?lang=pl) ma pierwszeństwo, potem język przeglądarki (Accept-Language: polski albo angielski),
# a gdy przeglądarka nie podaje żadnego z nich, używany jest LANG_DEFAULT; przełącza się przyciskiem PL | EN
const LANG_DEFAULT = "en"

# jednostki osi X: klucze (używane w kodzie) i napisy na liście
const UNIT_KEYS = ["nm", "Å", "µm", "cm⁻¹", "eV"]
const UNIT_OPTS_PL = ["nanometry (nm)", "angstremy (Å)", "mikrometry (µm)", "liczba falowa (cm⁻¹)", "energia fotonu (eV)"]
const UNIT_OPTS_EN = ["nanometers (nm)", "angstroms (Å)", "micrometers (µm)", "wavenumber (cm⁻¹)", "photon energy (eV)"]

"""
Przycina okno [lo, hi] [nm] do granic osi, zachowując jego szerokość (przy przesuwaniu zatrzymuje się
na brzegu). W skali logarytmicznej liczy w log10, więc zachowuje stosunek hi/lo. `positive = true`
(skala log albo energia/liczba falowa, gdzie λ = 0 nie istnieje) podnosi dolną granicę do X_MIN_LOG.
"""
function clamp_view(lo, hi, logmode; positive = logmode)
    t  = logmode ? log10 : identity
    ti = logmode ? (x -> 10.0^x) : identity
    if positive
        lo = max(lo, X_MIN_LOG)
        hi = max(hi, lo * 1.001)
    end
    amin = t(positive ? X_MIN_LOG : X_MIN)
    amax = t(X_MAX)
    a = t(lo)
    b = t(hi)
    w = clamp(b - a, logmode ? 0.002 : X_MIN_WIDTH, amax - amin)
    a = (a + b) / 2 - w / 2        # przy korekcie szerokości zachowaj środek okna (inaczej zoom ponad limit przesuwa widok)
    if a < amin
        a = amin
    end
    if a + w > amax
        a = amax - w
    end
    return (ti(a), ti(a + w))
end

# etykiety osi λ jako zwykłe liczby (100 zamiast 10², 1000 zamiast 10³)
function fmt_nm(vals)
    v = collect(Float64, vals)
    # liczba miejsc po przecinku z najmniejszego odstępu między tickami (żeby przy dużym zoomie się nie dublowały)
    d = length(v) >= 2 && minimum(abs.(diff(v))) > 0 ? clamp(ceil(Int, -log10(minimum(abs.(diff(v))))), 0, 8) : 3
    return [d == 0 ? string(round(Int, x)) : string(round(x; digits = d)) for x in v]
end

# "ładne" ticki osi λ dla okna [lo, hi] nm: 400, 500, 600... (liniowo lub 1-2-5 w skali log).
# Zwraca samą listę wartości (etykiety robi xtickformat = fmt_nm).
function nice_ticks(lo, hi, logmode)
    vals = Float64[]
    if logmode && hi / lo > 3
        mants = hi / lo > 10 ? (1, 2, 5) : (1, 2, 3, 5)
        for k in floor(Int, log10(lo)):ceil(Int, log10(hi)), m in mants
            v = m * 10.0^k
            lo <= v <= hi && push!(vals, v)
        end
    else
        raw  = (hi - lo) / 6
        mag  = 10.0^floor(log10(raw))
        step = mag * first(m for m in (1, 2, 5, 10) if m * mag >= raw)
        vals = collect(ceil(lo / step) * step : step : hi)
    end
    length(vals) < 2 && return Makie.automatic
    return vals
end

# PL: „1 linia spektralna”, „2 linie spektralne”, „5 linii spektralnych”; EN: „1 spectral line”, „2 spectral lines”
function lines_txt(n, lang)
    lang == "en" && return n == 1 ? "1 spectral line" : "$n spectral lines"
    n == 1 && return "1 linia spektralna"
    (2 <= n % 10 <= 4 && !(12 <= n % 100 <= 14)) && return "$n linie spektralne"
    return "$n linii spektralnych"
end

# względna suma natężeń: liczba całkowita od 10 wzwyż, poniżej 2 cyfry znaczące
fmt_rel(v) = !isfinite(v) ? "∞" : v >= 1e6 ? string(round(v; sigdigits = 3)) :
             v >= 10 ? string(round(Int, v)) : string(round(v; sigdigits = 2))

const ROMAN = ("I", "II", "III", "IV", "V")

# kolor kreski widoczny na jasnym tle: jasne barwy (żółty, cyjan) przyciemniamy do umiarkowanej jasności
function darken_for_light(c)
    lum = 0.299 * c[1] + 0.587 * c[2] + 0.114 * c[3]
    f = lum > 0.35 ? 0.35 / lum : 1.0           # próg dobrany tak, by kontrast na bieli wszędzie wynosił co najmniej ok. 3,8:1
    return (c[1] * f, c[2] * f, c[3] * f)
end

# udział procentowy: pełne procenty od 10% wzwyż, poniżej 2 cyfry znaczące
fmt_pct(f) = f >= 0.1 ? "$(round(Int, 100f))%" : "$(round(100f; sigdigits = 2))%"

# Odczyt ustawień z adresu strony (?z=26&T=8000&...): liczby i wartości logiczne z wartościami domyślnymi.
q_num(q, k, default) = (v = tryparse(Float64, get(q, k, "")); (v === nothing || !isfinite(v)) ? default : v)
q_int(q, k, default, lo, hi) = clamp(round(Int, q_num(q, k, default)), lo, hi)
q_bool(q, k, default) = haskey(q, k) ? q[k] == "1" : default

# Napisy interfejsu (klucz => (PL, EN)). Napisy zmieniające się dynamicznie (opis linii, pasek informacji,
# okno "O programie") są tłumaczone na miejscu przez L(pl, en) w build_app.
const TR = Dict(
    "element"   => ("Pierwiastek:", "Element:"),
    "symbol"    => ("Symbol:", "Symbol:"),
    "range"     => ("Zakres:", "Range:"),
    "view"      => ("Widok:", "View:"),
    "units"     => ("Jednostki:", "Units:"),
    "lines"     => ("Linie:", "Lines:"),
    "neutral"   => ("niezjonizowane", "neutral"),
    "ionized"   => ("zjonizowane:", "ionized:"),
    "settings"  => ("⚙ Ustawienia", "⚙ Settings"),
    "g_axis"    => ("Oś", "Axis"),
    "g_int"     => ("Natężenia", "Intensities"),
    "g_merge"   => ("Zlewanie linii", "Line merging"),
    "g_look"    => ("Wygląd", "Appearance"),
    "axscale"   => ("Skala osi:", "Axis scale:"),
    "medium"    => ("Długość fali:", "Wavelength:"),
    "source"    => ("Źródło:", "Source:"),
    "scale"     => ("Skala:", "Scale:"),
    "reference" => ("Odniesienie:", "Reference:"),
    "threshold" => ("Zakres jasności [dekady]:", "Brightness range [decades]:"),
    "thresholdtip" => ("Pokazuje linie słabsze od linii odniesienia co najwyżej o tyle rzędów wielkości (dekad); np. 4 = do 10 000 razy słabsze. " *
                       "W skali logarytmicznej ustala też skalę wysokości kresek.",
                       "Shows lines weaker than the reference line by at most this many orders of magnitude (decades); e.g. 4 = up to 10,000 times weaker. " *
                       "On the logarithmic scale it also sets the scale of the line heights."),
    "resol"     => ("Rozdzielczość [nm]:", "Resolution [nm]:"),
    "resolpx"   => ("rozdzielczość zależna od zoomu [px]:", "zoom-dependent resolution [px]:"),
    "saha"      => ("równowaga jonizacji (Saha)", "ionization equilibrium (Saha)"),
    "lg_ne"     => ("lg n_e [cm⁻³]:", "lg n_e [cm⁻³]:"),
    "bw"        => ("czarno-białe", "black and white"),
    "gray"      => ("linie poza światłem widzialnym szare", "lines outside visible light in gray"),
    "members"   => ("składowe zlanych linii w opisie", "components of merged lines in the description"),
    "tip1"      => ("atom niezjonizowany (X)", "neutral atom (X)"),
    "tip2"      => ("jednokrotnie zjonizowane (X⁺)", "singly ionized (X⁺)"),
    "tip3"      => ("dwukrotnie zjonizowane (X²⁺)", "doubly ionized (X²⁺)"),
    "tip4"      => ("trzykrotnie zjonizowane (X³⁺)", "triply ionized (X³⁺)"),
    "tip5"      => ("czterokrotnie zjonizowane (X⁴⁺)", "quadruply ionized (X⁴⁺)"),
    "langtip"   => ("Język interfejsu: PL / EN", "Interface language: PL / EN"),
    "save"      => ("💾 Zapisz obraz", "💾 Save image"),
    "presets"   => ("🧪 Przykłady", "🧪 Examples"),
    "presetnote" => ("Ustawienia poglądowe: pierwiastek, stopnie jonizacji, światło widzialne i temperatura dobrane tak, by pokazać charakterystyczne linie. " *
                     "W prawdziwych źródłach plazma często nie jest w równowadze LTE, więc proporcje natężeń mogą się różnić.",
                     "Illustrative settings: element, ionization stages, visible light and a temperature chosen to show the characteristic lines. " *
                     "In real sources the plasma is often not in LTE, so the intensity ratios may differ."),
    "linktip"   => ("Kopiuj link do bieżącego widoku (adres strony zawiera ustawienia, więc odświeżenie strony też je zachowuje)",
                    "Copy a link to the current view (the page address holds the settings, so reloading the page keeps them too)"),
    "tbl_title" => ("Linie w widoku", "Lines in view"),
    "warming"   => ("⏳ Rozgrzewanie…", "⏳ Warming up…"),
    "warmtip"   => ("Pierwsze użycie każdej funkcji kompiluje kod, co powoduje chwilowe lagi. Przycisk wykonuje przykładowe operacje " *
                    "(przesunięcia, powiększenia, skala, jednostki, suwaki temperatury i zakresu jasności, inny pierwiastek), potem przywraca ustawienia i znika.",
                    "The first use of each function compiles code, which causes brief lags. The button runs sample operations " *
                    "(panning, zooming, scale, units, the temperature and brightness range sliders, another element), then restores your settings and disappears."),
)

# Para list rozwijanych PL/EN o wspólnym indeksie (źródłem prawdy jest `idx`). Zmiana napisów na liście w
# Bonito resetuje wybór, więc zamiast tego trzymamy obie wersje i pokazujemy jedną z nich.
function dd2(opts_pl, opts_en; index = 1)
    a = Bonito.Dropdown(opts_pl; index = index)
    b = Bonito.Dropdown(opts_en; index = index)
    on(a.option_index) do i
        b.option_index[] == i || (b.option_index[] = i)
    end
    on(b.option_index) do i
        a.option_index[] == i || (a.option_index[] = i)
    end
    w = maximum(length, vcat(collect(String, opts_pl), collect(String, opts_en))) + 4     # szerokość w "ch"
    return (idx = a.option_index, a = a, b = b, w = w)
end

function build_app()
    return Bonito.App(title = "Spektruś.jl") do session, request
        q = parse_query(request.target)        # ustawienia z adresu strony (pusty słownik, gdy ich nie ma)
        t_load = time();  t_prev = Ref(t_load)
        load_lap(label) = LOAD_TIMING && (now = time();
            println("  [load] ", rpad(label, 36), lpad(round(Int, 1000 * (now - t_prev[])), 6), " ms   (łącznie ",
                    round(Int, 1000 * (now - t_load)), " ms)");
            t_prev[] = now; nothing)
        # ---------------- kontrolki ----------------
        # pierwiastek: nazwa, symbol i liczba atomowa to trzy powiązane kontrolki (jedna wspólna wartość Zsel)
        z0       = something(findfirst(==(q_int(q, "z", 11, 1, 99)), ELEMENTS.Z), findfirst(==(11), ELEMENTS.Z), 1)
        Zsel     = Observable(ELEMENTS.Z[z0])

        # język interfejsu i pomocnicze funkcje tłumaczenia
        accept = ""                              # nagłówek Accept-Language wysłany przez przeglądarkę
        for (k, v) in request.headers
            lowercase(String(k)) == "accept-language" && (accept = String(v); break)
        end
        lang0 = get(q, "lang", pick_language(accept; default = LANG_DEFAULT))
        lang = Observable(lang0 in ("pl", "en") ? lang0 : LANG_DEFAULT)
        L(pl, en) = lang[] == "en" ? en : pl                       # tekst w bieżącym języku
        tr(key) = map(l -> TR[key][l == "en" ? 2 : 1], lang)       # Observable napisu (zmienia się z językiem)
        show_for(l) = map(x -> x == l ? "display:flex;align-items:center;gap:6px;" : "display:none;", lang)
        # dwie wersje kontrolki (PL i EN), widoczna jedna; `w` = minimalna szerokość [ch], taka sama w obu językach
        pair(a, b; w = 0) = Bonito.DOM.span(Bonito.DOM.span(a; style = show_for("pl")),
                                            Bonito.DOM.span(b; style = show_for("en"));
                                            class = "pair",
                                            style = "display:inline-block;" * (w > 0 ? "min-width:$(round(0.5w; digits = 1))em;" : ""))
        ddv(d) = pair(d.a, d.b; w = d.w)
        # etykieta o stałej szerokości (najdłuższy z napisów PL/EN), żeby zmiana języka nie przesuwała kontrolek
        lbl(key) = Bonito.DOM.span(tr(key);
                                   style = "display:inline-block;min-width:$(round(0.48 * maximum(length, TR[key]); digits = 1))em;")

        name_d   = dd2(NAME_OPTS_PL, NAME_OPTS_EN; index = z0)
        sym_dd   = Bonito.Dropdown(SYM_OPTS; index = z0)
        # liczba atomowa: własne pole <input type=number> z Observable{Int} (Bonito.NumberInput istnieje tylko dla Float64)
        z_val    = Observable(Int(ELEMENTS.Z[z0]))
        z_in     = Bonito.DOM.input(; type = "number", value = z_val, min = 1, max = 99, step = 1,
                                    style = "width:3.5em;padding:2px 4px;font-weight:600;font-size:0.9em;border-radius:3px;",
                                    onchange = Bonito.js"""event => {
                                        const v = Math.round(parseFloat(event.srcElement.value));
                                        if (Number.isFinite(v) && $(z_val).value != v) { $(z_val).notify(v); }
                                    }""")
        unit_d   = dd2(UNIT_OPTS_PL, UNIT_OPTS_EN; index = q_int(q, "u", 1, 1, length(UNIT_OPTS_PL)))
        med_d    = dd2(["powietrze", "próżnia"], ["air", "vacuum"]; index = q_int(q, "med", 1, 1, 2))
        xs_d     = dd2(["liniowa", "logarytmiczna"], ["linear", "logarithmic"]; index = q_int(q, "xs", 1, 1, 2))
        scale_d  = dd2(["logarytmiczna", "liniowa"], ["logarithmic", "linear"]; index = q_int(q, "sc", 1, 1, 2))
        src_d    = dd2(["model LTE", "NIST"], ["LTE model", "NIST"]; index = q_int(q, "src", 1, 1, 2))
        ref_d    = dd2(["linia widzialna", "linia w oknie", "linia pierwiastka"],
                       ["visible-range line", "line in window", "element line"]; index = q_int(q, "ref", 1, 1, 3))
        ionset   = haskey(q, "ions") ? Set(filter(!isnothing, tryparse.(Int, split(q["ions"], ',')))) : Set([1])
        ion1_cb  = Bonito.Checkbox(1 in ionset)              # linie niezjonizowane (I)
        ion2_cb  = Bonito.Checkbox(2 in ionset)              # linie zjonizowane (II)
        ion3_cb  = Bonito.Checkbox(3 in ionset)              # podwójnie zjonizowane (III)
        ion4_cb  = Bonito.Checkbox(4 in ionset)              # potrójnie zjonizowane (IV)
        ion5_cb  = Bonito.Checkbox(5 in ionset)              # poczwórnie zjonizowane (V)
        saha_cb  = Bonito.Checkbox(q_bool(q, "saha", true))  # równowaga jonizacji (Saha)
        bw_cb    = Bonito.Checkbox(q_bool(q, "bw", false))   # tryb czarno-biały
        gray_cb  = Bonito.Checkbox(q_bool(q, "gray", false)) # linie UV i IR szare (jak w oryginalnym Spektrusiu)
        members_cb = Bonito.Checkbox(q_bool(q, "mem", false))    # lista składowych zlanych linii w opisie
        T_sl     = Bonito.Slider(1000:100:20000; value = snap(1000:100:20000, q_num(q, "T", 6000)))
        ne_sl    = Bonito.Slider(12.0:0.5:18.0; value = snap(12.0:0.5:18.0, q_num(q, "ne", 15.0)))     # lg n_e [cm⁻³] (dla Sahy)
        res_sl   = Bonito.Slider(0.0:0.1:5.0; value = snap(0.0:0.1:5.0, q_num(q, "res", 0.0)))
        thr_sl   = Bonito.Slider(0.5:0.1:6.0; value = snap(0.5:0.1:6.0, q_num(q, "thr", 4.0)))
        rescb    = Bonito.Checkbox(q_bool(q, "rescb", false))        # rozdzielczość zależna od zoomu
        px_sl    = Bonito.Slider(1:10; value = q_int(q, "px", 2, 1, 10))     # ... w pikselach ekranu
        b_left   = Bonito.Button("◀")
        b_right  = Bonito.Button("▶")
        b_zin    = Bonito.Button("＋")
        b_zout   = Bonito.Button("－")
        about_a  = Bonito.Button("ⓘ O programie")
        about_b  = Bonito.Button("ⓘ About")
        lang_btn = Bonito.Button("PL | EN")
        pngd_a   = Bonito.Button("Ciemne tło (jak na ekranie)")
        pngd_b   = Bonito.Button("Dark background (as on screen)")
        pngl_a   = Bonito.Button("Jasne tło (do druku)")
        pngl_b   = Bonito.Button("Light background (for printing)")
        warm_a   = Bonito.Button("🔥 Rozgrzej")
        warm_b   = Bonito.Button("🔥 Warm up")
        warm_state = Observable(WARMED[] ? "done" : "idle")        # "idle" -> "running" -> "done"
        about_open = Observable(false)
        table_a    = Bonito.Button("📋 Tabela linii")
        table_b    = Bonito.Button("📋 Line table")
        tbl_csv_a  = Bonito.Button("Pobierz CSV")
        tbl_csv_b  = Bonito.Button("Download CSV")
        tbl_tsv_a  = Bonito.Button("Pobierz TSV")
        tbl_tsv_b  = Bonito.Button("Download TSV")
        tbl_close_a = Bonito.Button("✕ Zamknij")
        tbl_close_b = Bonito.Button("✕ Close")
        link_btn   = Bonito.Button("🔗 Link")
        link_msg   = Observable("")
        table_open = Observable(false)
        table_rev  = Observable(0)           # zwiększany, gdy tabela ma się odświeżyć (zmiana widoku)
        help_a     = Bonito.Button("❓ Pomoc")
        help_b     = Bonito.Button("❓ Help")
        help_open  = Observable(false)
        # przykłady: (nazwa PL, nazwa EN, Z, T [K], stopnie jonizacji); T dobrane poglądowo (gazy o wysokim wzbudzeniu: wyższa)
        PRESETS = (("Lampa wodorowa", "Hydrogen lamp", 1, 10000, (1,)), ("Lampa helowa", "Helium lamp", 2, 10000, (1,)),
                   ("Lampa neonowa", "Neon lamp", 10, 10000, (1,)), ("Lampa sodowa", "Sodium lamp", 11, 6000, (1,)),
                   ("Lampa rtęciowa", "Mercury lamp", 80, 6000, (1,)),
                   ("Wapń: linie H i K Słońca (Ca II)", "Calcium: solar H and K lines (Ca II)", 20, 5800, (1, 2)))
        preset_btns = [(Bonito.Button(p[1]), Bonito.Button(p[2]), p[3], p[4], p[5]) for p in PRESETS]
        preset_busy = Ref(false)          # na czas stosowania przykładu wstrzymuje przerysowania (jedno na końcu)
        key_event  = Observable("")          # naciśnięty klawisz skrótu (przekazywany z przeglądarki, z numerem, by każde naciśnięcie było zmianą)

        info  = Observable("")
        note  = Observable("")            # ostrzeżenie pod paskiem informacji (np. stopień bez danych do LTE)
        hover = Observable(L("Najedź kursorem na linię.", "Hover over a line."))
        summary = Observable("")      # „N linii spektralnych o względnej sumie natężeń X”
        hint    = Observable("")      # podpowiedź o zaznaczaniu i powiększaniu
        members_txt = Observable("")  # składowe zlanych linii (gdy włączone)

        load_lap("kontrolki (widżety)")
        # ---------------- stan osi X (potrzebny już przy tworzeniu osi) ----------------
        log_x  = Ref(false)               # czy oś X jest logarytmiczna
        unit_x = Ref("nm")                # jednostka osi X (klucz z UNIT_KEYS)
        pos_now() = log_x[] || unit_is_energy(unit_x[])      # czy λ = 0 jest niedozwolone
        fmt_x(vals) = fmt_nm(log_x[] ? 10.0 .^ collect(Float64, vals) : vals)

        # ---------------- wykres ----------------
        fig = Figure(size = (1100, 520), backgroundcolor = RGBf(0.09, 0.09, 0.11))
        axstyle = (backgroundcolor = :black,
                   xlabelcolor = :white, ylabelcolor = :white,
                   xticklabelcolor = :white, yticklabelcolor = :white,
                   xtickcolor = :white, ytickcolor = :white,
                   leftspinecolor = :gray40, bottomspinecolor = :gray40,
                   topspinevisible = false, rightspinevisible = false,
                   xgridvisible = false, ygridvisible = false,
                   yzoomlock = true, ypanlock = true, yrectzoom = false)
        ax_strip = Axis(fig[1, 1]; height = 90,
                        yticksvisible = false, yticklabelsvisible = false,
                        xticksvisible = false, xticklabelsvisible = false, axstyle...)
        ax = Axis(fig[2, 1]; xlabel = "λ [nm]", ylabel = "natężenie względne",
                  xtickformat = fmt_x, axstyle...)
        linkxaxes!(ax_strip, ax)
        xlims!(ax, 380, 780)
        load_lap("wykres (Figure, osie)")

        # ---------------- stan ----------------
        xq        = tryparse.(Float64, split(get(q, "x", ""), ','))
        x_init    = (length(xq) == 2 && all(!isnothing, xq) && all(isfinite, xq) && xq[1] < xq[2]) ? (xq[1], xq[2]) : (380.0, 780.0)
        view_x    = Ref(x_init)           # bieżący zakres osi X [nm] (w nm, nie w log10); zacznij od widoku z adresu
        redrawing = Ref(false)            # ignorujemy zmiany limitów wywołane przez nasz kod
        print_mode = Ref(false)           # tryb zapisu na jasnym tle: tylko na czas zrzutu PNG
        cache     = Ref{Any}(nothing)     # ostatnio policzone linie (element, T, ośrodek, źródło)
        shown     = Ref{Any}(nothing)     # to, co widać na wykresie (dla podpowiedzi)
        cur_res   = Ref(RES_FLOOR)        # rozdzielczość zlewania użyta przy ostatnim rysowaniu [nm]

        # ---------------- adres strony = stan widoku (zapamiętywanie ustawień i link) ----------------
        # Adres w przeglądarce jest aktualizowany po każdej zmianie widoku (z opóźnieniem 0,6 s), więc odświeżenie
        # strony przywraca ustawienia, a skopiowany adres otwiera dokładnie ten sam widok.
        function current_qs()
            x0, x1 = view_x[]
            ions = [z for (z, cb) in enumerate((ion1_cb, ion2_cb, ion3_cb, ion4_cb, ion5_cb)) if cb.value[]]
            kv = ("z" => Zsel[], "x" => "$(round(x0; sigdigits = 8)),$(round(x1; sigdigits = 8))",
                  "u" => unit_d.idx[], "xs" => xs_d.idx[], "sc" => scale_d.idx[], "src" => src_d.idx[],
                  "ref" => ref_d.idx[], "med" => med_d.idx[], "ions" => join(ions, ','),
                  "T" => T_sl.value[], "ne" => ne_sl.value[], "thr" => thr_sl.value[], "res" => res_sl.value[],
                  "px" => px_sl.value[], "rescb" => Int(rescb.value[]), "saha" => Int(saha_cb.value[]),
                  "bw" => Int(bw_cb.value[]), "gray" => Int(gray_cb.value[]), "mem" => Int(members_cb.value[]), "lang" => lang[])
            return join(["$k=$v" for (k, v) in kv], "&")
        end

        # tabela pojedynczych linii w bieżącym widoku (z progiem jasności jak na wykresie), od najsilniejszej
        function current_table()
            s = shown[]
            s === nothing && return nothing
            x0, x1 = view_x[]
            rows = Int[]
            for i in eachindex(s.first)
                append!(rows, s.first[i]:s.last[i])
            end
            rows = [r for r in rows if x0 <= s.lines.λ[r] <= x1]
            isempty(rows) && return nothing
            return line_table(s.lines, rows, s.Iref, unit_x[])
        end

        url_ready  = Ref(false)           # adres aktualizujemy dopiero po pierwszym rysowaniu
        url_last   = Ref("")
        sync_wait  = Ref(false)
        function ui_sync!()               # zlewa serię zmian widoku w jedną aktualizację adresu i tabeli
            sync_wait[] && return nothing
            sync_wait[] = true
            @async begin
                sleep(0.6)
                sync_wait[] = false
                try
                    if url_ready[]
                        qs = current_qs()
                        if qs != url_last[]
                            url_last[] = qs
                            Bonito.evaljs(session, Bonito.js"""history.replaceState(null, '', '?' + $(qs));""")
                        end
                    end
                    table_open[] && (table_rev[] += 1)
                catch e
                    @warn "Aktualizacja adresu/tabeli nie powiodła się" exception = (e, catch_backtrace())
                end
            end
            return nothing
        end

        # --- współrzędne osi X ---
        # Oś Makie jest zawsze liniowa; skalę logarytmiczną i jednostki realizujemy sami: X to wartość w
        # jednostce `unit_x[]` (albo jej log10 przy skali log). Dzięki temu limity, pozycja myszy i ticki
        # są jednoznaczne (bez przeliczeń po stronie Makie). Okno widoku `view_x` zostaje w nm.
        to_X(λ) = (u = unit_value(λ, unit_x[]); log_x[] ? log10(u) : u)
        function from_X_nm(X)
            u = log_x[] ? 10.0^X : X
            unit_is_energy(unit_x[]) && (u = max(u, 1e-12))     # energia <= 0 oznacza λ -> ∞; przytnie to clamp_view
            return unit_to_nm(u, unit_x[])
        end
        function view_to_X(v)                       # okno w nm -> (X_min, X_max)
            a, b = to_X(v[1]), to_X(v[2])
            return a <= b ? (a, b) : (b, a)
        end
        function X_to_view(a, b)                    # (X, X) -> okno w nm (lo, hi)
            n1, n2 = from_X_nm(a), from_X_nm(b)
            return (min(n1, n2), max(n1, n2))
        end
        function ticks_for(lo, hi)                  # ticki w jednostce osi, jako pozycje X
            ua, ub = unit_value(lo, unit_x[]), unit_value(hi, unit_x[])
            vals = nice_ticks(min(ua, ub), max(ua, ub), log_x[])
            vals === Makie.automatic && return vals
            return log_x[] ? log10.(vals) : vals
        end

        last_ticks = Ref{Any}(nothing)
        function set_ticks!(t)                       # nie przebudowuj tick-ów, gdy się nie zmieniły
            t == last_ticks[] && return nothing
            last_ticks[] = t
            ax.xticks = t
            return nothing
        end

        # przerysowanie po zakończeniu zoomu/przesuwania (czekamy 0,2 s na brak kolejnych zmian)
        last_change = Ref(0.0)
        function schedule_redraw!()
            stamp = time()
            last_change[] = stamp
            @async begin
                sleep(0.2)
                last_change[] == stamp && request_redraw!()
            end
            return nothing
        end

        # Kreski tworzymy raz i potem tylko podmieniamy dane (tworzenie i usuwanie wykresów w WGLMakie
        # jest drogie). Liczba wierzchołków to stała pojemność `cap`; nadmiar to przezroczyste odcinki
        # zerowej długości w punkcie (0, 0). Dzięki temu tablice punktów
        # i kolorów mają zawsze tę samą długość. Wykresy budujemy od nowa tylko wtedy, gdy zabraknie pojemności.
        segs  = Ref{Any}(nothing)
        PAD   = Point2d(0.0, 0.0)
        BLANK = RGBAf(0, 0, 0, 0)

        # nakładka pokazująca zaznaczany obszar (przeciąganie lewym przyciskiem); x w przestrzeni osi X
        sel       = Observable((0.0, 0.0))
        sel_vis   = Observable(false)
        sel_main  = lift(sel) do r
            Rect2d(Point2d(min(r[1], r[2]), 0), Vec2d(abs(r[2] - r[1]), 1.05))
        end
        sel_strip = lift(sel) do r
            Rect2d(Point2d(min(r[1], r[2]), 0), Vec2d(abs(r[2] - r[1]), 1))
        end

        function drop_segments!()
            empty!(ax)
            empty!(ax_strip)
            segs[] = nothing
            return nothing
        end

        function make_segments!(cap)
            pts1 = Observable(fill(PAD, 2 * cap));  col1 = Observable(fill(BLANK, 2 * cap))
            pts2 = Observable(fill(PAD, 2 * cap));  col2 = Observable(fill(BLANK, 2 * cap))
            linesegments!(ax_strip, pts1; color = col1, linewidth = 1.5)
            linesegments!(ax, pts2; color = col2, linewidth = 1.5)
            for (axis, rect) in ((ax_strip, sel_strip), (ax, sel_main))
                poly!(axis, rect; color = RGBAf(1, 1, 1, 0.15), strokecolor = RGBAf(1, 1, 1, 0.9),
                      strokewidth = 1.5, visible = sel_vis, inspectable = false)
            end
            segs[] = (; pts1, col1, pts2, col2, cap)
            return nothing
        end

        # pomiar czasu etapów redraw! (wypisywany w konsoli; ustaw TIMING = false, żeby wyłączyć)
        TIMING = false
        tmark = Ref(time())
        function lap(label)
            TIMING || return nothing
            t = time()
            println("  [redraw] ", rpad(label, 34), lpad(round(Int, 1000 * (t - tmark[])), 6), " ms")
            tmark[] = t
            return nothing
        end

        # przypisuj podpowiedź tylko, gdy tekst się zmienił (Observable rozsyła każde przypisanie)
        set_hover!(txt) = (hover[] == txt || (hover[] = txt); nothing)
        set_members!(txt) = (members_txt[] == txt || (members_txt[] = txt); nothing)
        # podpowiedź pod opisem zlanych linii: tylko wtedy, gdy coś jeszcze da się rozdzielić
        function merge_hint(n, spread)
            n > 1 || return ""
            spread == 0 && return L("Te linie mają identyczną długość fali w danych, więc nie da się ich rozdzielić.",
                                    "These lines have identical wavelengths in the data, so they cannot be resolved.")
            cur_res[] <= RES_FLOOR * 1.001 && return ""          # najdrobniejsza rozdzielczość: nic już nie rozdzieli
            rescb.value[] && return L("Możesz zaznaczyć i powiększyć ten wycinek w celu zwiększenia rozdzielczości.",
                                      "You can select and zoom into this section to increase the resolution.")
            return L("Zmniejsz rozdzielczość (suwak) albo włącz jej zależność od zoomu, żeby rozdzielić te linie.",
                     "Decrease the resolution (slider) or make it zoom-dependent to resolve these lines.")
        end
        function set_summary!(sum_txt, hint_txt)
            summary[] == sum_txt || (summary[] = sum_txt)
            hint[] == hint_txt || (hint[] = hint_txt)
            return nothing
        end

        function redraw!()
            t_start = time(); tmark[] = t_start
            redrawing[] = true
            try
                Z      = Zsel[]
                medium = med_d.idx[] == 1 ? "air" : "vac"
                source = src_d.idx[] == 1 ? "lte" : "nist"
                scale  = scale_d.idx[] == 1 ? "log" : "lin"
                refmode = ref_d.idx[] == 2 ? "window" :
                          ref_d.idx[] == 3 ? "all" : "visible"
                T       = Float64(T_sl.value[])
                decades = Float64(thr_sl.value[])
                ions    = Int[]
                ion1_cb.value[] && push!(ions, 1)
                ion2_cb.value[] && push!(ions, 2)
                ion3_cb.value[] && push!(ions, 3)
                ion4_cb.value[] && push!(ions, 4)
                ion5_cb.value[] && push!(ions, 5)
                saha    = saha_cb.value[]
                ne_exp  = Float64(ne_sl.value[])
                bw      = bw_cb.value[]
                gray_out = gray_cb.value[]
                logmode = xs_d.idx[] == 2
                unit_key = UNIT_KEYS[clamp(unit_d.idx[], 1, length(UNIT_KEYS))]
                energy   = unit_is_energy(unit_key)
                energy && (medium = "vac")      # cm⁻¹ i eV: zawsze długości fali w próżni
                log_x[]  = logmode
                unit_x[] = unit_key

                key = (Z, Tuple(ions), medium, source, T, saha, saha ? ne_exp : 0.0)
                if cache[] === nothing || cache[].key != key
                    spec = prepare_lines(DF, Z; ions = ions, medium = medium, source = source, T = T,
                                         saha = saha, ne = 10.0^ne_exp)
                    cache[] = (key = key, spec = spec)
                end
                lines   = cache[].spec.lines
                use_lte = cache[].spec.use_lte
                skipped = cache[].spec.skipped
                lap("odczyt ustawień + linie (cache)")

                x0, x1 = clamp_view(view_x[]..., logmode; positive = logmode || energy)
                view_x[] = (x0, x1)

                res = max(Float64(res_sl.value[]), RES_FLOOR)
                if rescb.value[]
                    # rozdzielczość = px_sl pikseli ekranu przy bieżącym powiększeniu
                    wpx = max(Float64(ax.scene.viewport[].widths[1]), 200.0)
                    px  = Float64(px_sl.value[])
                    res = max(logmode ? px * log(x1 / x0) / wpx * sqrt(x0 * x1) :
                                        px * (x1 - x0) / wpx, RES_FLOOR)
                end
                cur_res[] = res
                cl = cluster_window(lines, -Inf, Inf, res)
                lap("cluster_window ($(length(lines.I)) linii)")

                Iref = reference_intensity(lines, x0, x1; ref = refmode)
                idx  = isempty(cl.I) ? Int[] : findall(visible_mask(cl.I, Iref, decades))

                nseg = length(idx)
                if segs[] === nothing || nseg > segs[].cap
                    drop_segments!()
                    make_segments!(max(2 * nseg, 4096))
                    lap("nowe wykresy (pojemność $(segs[].cap))")
                end
                sg = segs[]
                pts1 = fill(PAD, 2 * sg.cap);  col1 = fill(BLANK, 2 * sg.cap)
                pts2 = fill(PAD, 2 * sg.cap);  col2 = fill(BLANK, 2 * sg.cap)
                for (j, k) in enumerate(idx)
                    h = line_height(cl.I[k], Iref, decades; scale = scale)
                    X = to_X(cl.λ[k])
                    (cr, cg, cb) = bw ? (1.0, 1.0, 1.0) :
                                   (gray_out && !(380.0 <= cl.λ[k] <= 780.0)) ? (0.62, 0.62, 0.62) : wl2rgb(cl.λ[k])
                    # zapis na jasnym tle (do druku): czarne albo przyciemnione kreski, w głównym wykresie i na pasku
                    (mr, mg, mb) = !print_mode[] ? (cr, cg, cb) : bw ? (0.0, 0.0, 0.0) : darken_for_light((cr, cg, cb))
                    a = print_mode[] ? 0.3 + 0.7 * h : 0.15 + 0.85 * h        # na bieli słabe linie bledną, więc podnosimy minimum
                    i2 = 2 * j
                    pts1[i2 - 1] = Point2d(X, 0.0);  pts1[i2] = Point2d(X, 1.0)
                    col1[i2 - 1] = RGBAf(mr, mg, mb, a);  col1[i2] = RGBAf(mr, mg, mb, a)
                    pts2[i2 - 1] = Point2d(X, 0.0);  pts2[i2] = Point2d(X, h)
                    col2[i2 - 1] = RGBAf(mr, mg, mb, 1);  col2[i2] = RGBAf(mr, mg, mb, 1)
                end
                lap("punkty i kolory ($nseg kresek)")
                sg.pts1[] = pts1;  sg.col1[] = col1
                sg.pts2[] = pts2;  sg.col2[] = col2
                lap("aktualizacja danych wykresów")
                xlabel_txt = unit_key == "cm⁻¹" ? L("liczba falowa σ = 1/λ [cm⁻¹]", "wavenumber σ = 1/λ [cm⁻¹]") :
                             unit_key == "eV" ? L("energia fotonu [eV]", "photon energy [eV]") :
                             "λ [$unit_key, $(medium == "air" ? L("powietrze", "air") : L("próżnia", "vacuum"))]"
                ax.xlabel[] == xlabel_txt || (ax.xlabel = xlabel_txt)
                ylabel_txt = L("natężenie względne", "relative intensity")
                ax.ylabel[] == ylabel_txt || (ax.ylabel = ylabel_txt)
                xa, xb = view_to_X((x0, x1))
                xlims!(ax, xa, xb)
                set_ticks!(ticks_for(x0, x1))
                ylims!(ax, 0, 1.05)
                ylims!(ax_strip, 0, 1)
                lap("limity i ticki")

                shown[] = (λ = cl.λ[idx], I = cl.I[idx], Iref = Iref, rows = cl.rows[idx], n = cl.n[idx], spread = cl.spread[idx],
                           first = cl.first[idx], last = cl.last[idx], lines = lines)
                frac = cache[].spec.fractions
                saha_txt = frac === nothing ? "" :
                    ", Saha: lg n_e = $(round(ne_exp; digits = 1)), $(L("udziały", "fractions")) " *
                    join(["$(ROMAN[z]) $(fmt_pct(frac[z]))" for z in eachindex(frac) if frac[z] >= 0.0005], ", ")
                src_txt = isempty(lines.I) ? L("brak linii (dla wybranych stopni jonizacji nie ma danych)",
                                               "no lines (no data for the selected ionization stages)") :
                          use_lte ? L("model LTE", "LTE model") * ", T = $(Int(T)) K$saha_txt" :
                          (source == "lte" ? L("NIST (brak Aki dla tego pierwiastka, LTE niedostępny)",
                                               "NIST (no Aki for this element, LTE unavailable)") : "NIST")
                note[] = isempty(skipped) ? "" :
                    L("⚠ Model LTE nie narysuje linii stopnia " * join(["$(ROMAN[z]) ($(skipped[z]) linii)" for z in sort!(collect(keys(skipped)))], ", ") *
                      ": NIST nie podaje dla nich współczynnika Aki. Wybierz źródło „NIST”, żeby je zobaczyć " *
                      "(natężenia umowne, niewspółmierne między stopniami).",
                      "⚠ The LTE model will not draw the lines of stage " * join(["$(ROMAN[z]) ($(skipped[z]) lines)" for z in sort!(collect(keys(skipped)))], ", ") *
                      ": NIST gives no Aki coefficient for them. Choose the “NIST” source to see them " *
                      "(arbitrary intensities, not comparable between stages).")
                info[] = L("Linie widoczne: $(length(idx)) z $(length(cl.λ)) (całe widmo)   ·   " *
                           "źródło natężeń: $src_txt   ·   odniesienie jasności: $(round(Iref; sigdigits = 3))   ·   " *
                           "rozdzielczość: $(round(res; sigdigits = 2)) nm",
                           "Visible lines: $(length(idx)) of $(length(cl.λ)) (whole spectrum)   ·   " *
                           "intensity source: $src_txt   ·   brightness reference: $(round(Iref; sigdigits = 3))   ·   " *
                           "resolution: $(round(res; sigdigits = 2)) nm")
                lap("info")
                ui_sync!()
                TIMING && println("  [redraw] ", rpad("RAZEM (strona Julii)", 34),
                                  lpad(round(Int, 1000 * (time() - t_start)), 6), " ms")
            finally
                redrawing[] = false
            end
            return nothing
        end

        # Zlewanie przerysowań: gdy trwa rysowanie, kolejne zmiany (np. przeciąganie suwaka T, który wysyła
        # zdarzenie przy każdym kroku) tylko oznaczają stan jako "brudny"; po zakończeniu rysujemy jeszcze raz,
        # już dla ostatniego stanu. Dzięki temu nie powstaje kolejka przerysowań, którą trzeba by opróżniać.
        busy  = Threads.Atomic{Bool}(false)
        dirty = Threads.Atomic{Bool}(false)
        function request_redraw!()
            dirty[] = true
            Threads.atomic_cas!(busy, false, true) == false || return nothing   # ktoś już rysuje
            try
                while dirty[]
                    dirty[] = false
                    redraw!()
                end
            finally
                busy[] = false
            end
            dirty[] && request_redraw!()        # zmiana, która wpadła tuż przed zwolnieniem blokady
            return nothing
        end

        # ---------------- reakcje ----------------
        onany(Zsel, unit_d.idx, med_d.idx, xs_d.idx, scale_d.idx, src_d.idx, ref_d.idx,
              ion1_cb.value, ion2_cb.value, ion3_cb.value, ion4_cb.value, ion5_cb.value,
              saha_cb.value, ne_sl.value, bw_cb.value, gray_cb.value, T_sl.value, res_sl.value, thr_sl.value, rescb.value, px_sl.value) do args...
            preset_busy[] && return nothing
            request_redraw!()
        end

        # powiązane kontrolki pierwiastka: nazwa, symbol i liczba atomowa zmieniają się razem
        # (jedna wspólna wartość Zsel; sprawdzanie równości zapobiega pętlom)
        function set_element!(Z)
            Z == Zsel[] || (Zsel[] = Z)
            return nothing
        end
        on(name_d.idx) do i
            set_element!(ELEMENTS.Z[i])
        end
        on(sym_dd.option_index) do i
            set_element!(ELEMENTS.Z[i])
        end
        on(z_val) do v
            (v isa Integer || (v isa Real && isfinite(v) && isinteger(v))) && Int(v) in ELEMENTS.Z && set_element!(Int(v))
        end
        on(Zsel) do Z
            i = findfirst(==(Z), ELEMENTS.Z)
            i === nothing && return
            name_d.idx[] == i || (name_d.idx[] = i)
            sym_dd.option_index[] == i || (sym_dd.option_index[] = i)
            z_val[] == Z || (z_val[] = Int(Z))
        end

        # powiększanie/pomniejszanie (factor) i przesuwanie (shift, ułamek nowej szerokości okna);
        # liczymy w przestrzeni osi X, więc ▶ zawsze przesuwa w prawo na ekranie (także w cm⁻¹ i eV)
        function move_view!(factor, shift)
            xa, xb = view_to_X(view_x[])
            w = (xb - xa) * factor
            a = (xa + xb) / 2 - w / 2 + shift * w
            lo, hi = X_to_view(a, a + w)
            view_x[] = clamp_view(lo, hi, log_x[]; positive = pos_now())
            request_redraw!()
        end
        on(b_zin.value) do _
            move_view!(0.5, 0.0)
        end
        on(b_zout.value) do _
            move_view!(2.0, 0.0)
        end
        on(b_left.value) do _
            move_view!(1.0, -0.25)
        end
        on(b_right.value) do _
            move_view!(1.0, 0.25)
        end

        # okno „O programie” (przełączane przyciskiem)
        for b in (about_a, about_b)
            on(b.value) do _
                about_open[] = !about_open[]
                about_open[] && (help_open[] = false)
            end
        end
        for b in (help_a, help_b)
            on(b.value) do _
                help_open[] = !help_open[]
                help_open[] && (about_open[] = false)
            end
        end
        on(lang_btn.value) do _
            lang[] = lang[] == "en" ? "pl" : "en"
        end
        # opis funkcji podziału zależy od tego, czy jest levels.csv (pełne listy poziomów NIST)
        levels_full = level_db() !== nothing
        U_pl = levels_full ?
            "Funkcje podziału U_z liczone są z pełnych list poziomów NIST (levels.csv), uzupełnionych o poziomy z pobranych linii; poziomy poniżej energii jonizacji. " :
            "Funkcje podziału U_z liczone są z poziomów występujących w pobranych liniach, więc to przybliżenie (dla złożonych atomów może zaniżać U). "
        U_en = levels_full ?
            "The partition functions U_z are computed from the full NIST level lists (levels.csv), supplemented with levels from the downloaded lines; levels below the ionization energy. " :
            "The partition functions U_z are computed from the levels that appear in the downloaded lines, so this is an approximation (for complex atoms it may underestimate U). "
        about_panel = map(about_open, lang) do o, _l
            o || return Bonito.DOM.div()
            Bonito.DOM.div(
                Bonito.DOM.h3("Spektruś.jl"; style = "margin:0 0 6px 0;"),
                Bonito.DOM.p(Bonito.DOM.b(L("Autorzy tej wersji: ", "Authors of this version: ")), AUTHOR,
                             L(" oraz Claude (Anthropic).", " and Claude (Anthropic).")),
                Bonito.DOM.p(L("Kolory linii: przybliżenie długości fali na RGB według algorytmu Dana Brutona (tak jak w oryginalnym Spektrusiu); " *
                               "poza światłem widzialnym kolory są umowne.",
                               "Line colors: wavelength-to-RGB approximation by Dan Bruton's algorithm (as in the original Spektruś); " *
                               "outside visible light the colors are conventional.")),
                Bonito.DOM.p(L("Kod źródłowy: ", "Source code: "),
                             Bonito.DOM.a(REPO_URL; href = REPO_URL, target = "_blank", style = "color:#8cf;")),
                Bonito.DOM.p(L("Odtworzenie programu Spektruś 1.0 (Mikołaj Pytel, 2002) w Julii. " *
                               "Pokazuje widma emisyjne pierwiastków 1–99 w UV, świetle widzialnym i podczerwieni.",
                               "A recreation of the program Spektruś 1.0 (Mikołaj Pytel, 2002) in Julia. " *
                               "It shows the emission spectra of elements 1–99 in the UV, visible and infrared.")),
                Bonito.DOM.p(L("Dane linii spektralnych: ", "Spectral line data: "), Bonito.DOM.b("NIST Atomic Spectra Database"),
                             L(" (stopnie jonizacji I–V, 100–2500 nm). Prosimy cytować: A. Kramida, Yu. Ralchenko, J. Reader " *
                               "i NIST ASD Team, NIST Atomic Spectra Database, ",
                               " (ionization stages I–V, 100–2500 nm). Please cite: A. Kramida, Yu. Ralchenko, J. Reader " *
                               "and NIST ASD Team, NIST Atomic Spectra Database, "),
                             Bonito.DOM.a("https://physics.nist.gov/asd"; href = "https://physics.nist.gov/asd",
                                          target = "_blank", style = "color:#8cf;"), "."),
                Bonito.DOM.p(L("Natężenia: model LTE, emisja optycznie cienka, I ∝ (N_z/N)·g_k·A_ki·ν·exp(−E_k/kT)/U_z(T) " *
                               "(temperatura z suwaka T). Udziały stopni jonizacji N_z/N wynikają z równowagi Saha dla " *
                               "gęstości elektronowej n_e (suwak lg n_e); uwzględniono stopnie I–V. " * U_pl * "Energie jonizacji: tabela z pakietu mendeleev (zgodna z NIST ASD w sprawdzonych " *
                               "przypadkach). Po wyłączeniu opcji Saha jasności linii różnych stopni jonizacji nie są " *
                               "porównywalne, a w obrębie jednego stopnia są.",
                               "Intensities: LTE model, optically thin emission, I ∝ (N_z/N)·g_k·A_ki·ν·exp(−E_k/kT)/U_z(T) " *
                               "(temperature from the T slider). The ionization stage fractions N_z/N come from the Saha equilibrium " *
                               "for the electron density n_e (lg n_e slider); stages I–V are included. " * U_en * "Ionization energies: a table from the mendeleev package " *
                               "(consistent with NIST ASD in the cases checked). With the Saha option off, the brightnesses of lines " *
                               "of different ionization stages are not comparable, while within a single stage they are.")),
                Bonito.DOM.p(Bonito.DOM.b(L("Źródła modelu fizycznego. ", "Sources of the physical model. ")),
                             L("Rozkład Boltzmanna populacji poziomów, emisja optycznie cienka w LTE oraz równowaga jonizacji Saha " *
                               "to standardowe elementy fizyki plazmy; opisano je m.in. w:",
                               "The Boltzmann distribution of level populations, optically thin emission in LTE and the Saha " *
                               "ionization equilibrium are standard elements of plasma physics; they are described, among others, in:")),
                Bonito.DOM.ul(
                    Bonito.DOM.li("H. R. Griem, ", Bonito.DOM.i("Principles of Plasma Spectroscopy"),
                                  L(", Cambridge University Press, 1997 (rozdz. 7: Thermodynamic equilibrium relations), ",
                                    ", Cambridge University Press, 1997 (ch. 7: Thermodynamic equilibrium relations), "),
                                  Bonito.DOM.a(L("strona książki", "book page");
                                               href = "https://www.cambridge.org/core/books/principles-of-plasma-spectroscopy/90F6DBB512089B7E3AC121744EFE3D93",
                                               target = "_blank", style = "color:#8cf;"), "."),
                    Bonito.DOM.li("M. N. Saha, ", Bonito.DOM.i("Ionization in the solar chromosphere"),
                                  ", Philosophical Magazine 40, 472 (1920), doi: ",
                                  Bonito.DOM.a("10.1080/14786441008636148";
                                               href = "https://doi.org/10.1080/14786441008636148",
                                               target = "_blank", style = "color:#8cf;"), ".");
                    style = "margin:4px 0 8px 0;padding-left:20px;"),
                Bonito.DOM.p(L("Stała w równaniu Saha: 2(2π m_e k / h²)^{3/2} = 4,8294·10¹⁵ cm⁻³ K⁻³ᐟ² (stałe CODATA). " *
                               "Model jest orientacyjny: zakłada LTE i optycznie cienką plazmę, które nie zawsze są spełnione.",
                               "Constant in the Saha equation: 2(2π m_e k / h²)^{3/2} = 4.8294·10¹⁵ cm⁻³ K⁻³ᐟ² (CODATA constants). " *
                               "The model is approximate: it assumes LTE and an optically thin plasma, which are not always satisfied.")),
                Bonito.DOM.p(L("Obsługa: kółko myszy = zoom, przeciągnięcie lewym = zaznaczenie obszaru, " *
                               "prawym = przesuwanie; przyciski ◀ ▶ ＋ － działają bez myszy. Jednostki osi: nm, Å, µm, cm⁻¹, eV " *
                               "(cm⁻¹ i eV nie zależą od ośrodka: liczone są z długości fali w próżni, a oś rośnie z energią; " *
                               "przełącznik powietrze/próżnia jest wtedy ukryty).",
                               "Controls: mouse wheel = zoom, drag with the left button = select an area to zoom into, " *
                               "right button = pan; the buttons ◀ ▶ ＋ － work without a mouse. Axis units: nm, Å, µm, cm⁻¹, eV " *
                               "(cm⁻¹ and eV do not depend on the medium: they are computed from vacuum wavelengths and the axis " *
                               "increases with energy; the air/vacuum switch is then hidden)."));
                style = "max-width:760px;margin:8px 0;padding:10px 14px;background:#222228;" *
                        "border:1px solid #444;border-radius:6px;font-size:0.92em;line-height:1.4;" *
                        "box-sizing:border-box;max-height:70vh;overflow:auto;")
        end

        range_nodes = Any[]
        w_range = maximum(max(length(r[1]), length(r[2])) for r in RANGES) + 3     # jedna szerokość dla wszystkich czterech
        for (name_pl, name_en, (lo, hi)) in RANGES
            ba = Bonito.Button(name_pl)
            bb = Bonito.Button(name_en)
            for b in (ba, bb)
                on(b.value) do _
                    view_x[] = (lo, hi)
                    request_redraw!()
                end
            end
            push!(range_nodes, pair(ba, bb; w = w_range))
        end

        # zakres ustawiony myszą (zoom/przesuwanie): zapamiętujemy go i pilnujemy granic osi
        on(ax.finallimits) do r
            redrawing[] && return
            lo, hi = X_to_view(r.origin[1], r.origin[1] + r.widths[1])
            clo, chi = clamp_view(lo, hi, log_x[]; positive = pos_now())
            if abs(clo - lo) > 1e-9 * max(1.0, abs(lo)) || abs(chi - hi) > 1e-9 * max(1.0, abs(hi))
                redrawing[] = true
                try
                    xa, xb = view_to_X((clo, chi))
                    xlims!(ax, xa, xb)
                finally
                    redrawing[] = false
                end
            end
            view_x[] = (clo, chi)
            set_ticks!(ticks_for(clo, chi))
            ui_sync!()
            (ref_d.idx[] == 2 || rescb.value[]) && schedule_redraw!()
        end

        # podpowiedź po najechaniu na linię
        cursor = Observable("")                       # odczyt położenia kursora w jednostce osi (niezależnie od linii)
        cursor_digits = Dict("nm" => 3, "Å" => 2, "µm" => 5, "cm⁻¹" => 1, "eV" => 4)
        function hover_handler(_)
            s = shown[]
            s === nothing && return
            if !(Makie.is_mouseinside(ax.scene) || Makie.is_mouseinside(ax_strip.scene))
                cursor[] = ""
                return
            end
            let u = unit_x[], v = unit_value(from_X_nm(mouseposition(ax.scene)[1]), unit_x[])
                cursor[] = L("kursor: ", "cursor: ") * "$(round(v; digits = get(cursor_digits, u, 3))) $u"
            end
            λs = s.λ
            if isempty(λs)
                set_hover!(L("Brak widocznych linii.", "No visible lines."))
                set_summary!("", "")
                set_members!("")
                return
            end
            mX = mouseposition(ax.scene)[1]       # oś Makie jest zawsze liniowa: to współrzędna X
            mx = from_X_nm(mX)                    # to samo w nm
            xa, xb = view_to_X(view_x[])
            tol = 6 * (xb - xa) / ax.scene.viewport[].widths[1]      # ok. 6 pikseli
            k = searchsortedfirst(λs, mx)         # λs rośnie z nm niezależnie od jednostki osi
            cand = [i for i in (k - 1, k) if 1 <= i <= length(λs)]
            i = cand[argmin([abs(to_X(λs[j]) - mX) for j in cand])]
            if abs(to_X(λs[i]) - mX) <= tol
                set_hover!(describe_line(s.lines, s.rows[i]; n = s.n[i], spread = s.spread[i], unit = unit_x[], lang = lang[]))
                set_summary!("$(lines_txt(s.n[i], lang[])) $(L("o względnej sumie natężeń", "with a relative intensity sum of")) $(fmt_rel(100 * s.I[i] / s.Iref))",
                             merge_hint(s.n[i], s.spread[i]))
                set_members!(members_cb.value[] && s.n[i] > 1 ? describe_members(s.lines, s.first[i], s.last[i]; lang = lang[]) : "")
            else
                set_hover!(L("Najedź kursorem na linię.", "Hover over a line."))
                set_summary!("", "")
                set_members!("")
            end
        end
        on(events(fig).mouseposition) do ev
            try
                hover_handler(ev)
            catch e
                @warn "Błąd w podpowiedzi (hover)" exception = (e, catch_backtrace())
            end
        end

        # podgląd zaznaczenia: rysujemy prostokąt od miejsca naciśnięcia lewego przycisku do kursora
        # (samo powiększenie po puszczeniu przycisku robi Makie)
        dragging = Ref(false)
        drag_x0  = Ref(0.0)
        mouse_X() = mouseposition(ax.scene)[1]      # oś Makie jest zawsze liniowa: współrzędna X
        in_axes() = Makie.is_mouseinside(ax.scene) || Makie.is_mouseinside(ax_strip.scene)
        # wysoki priorytet: Makie przechwytuje zdarzenia myszy w trakcie przeciągania (własny zoom),
        # więc listenery o priorytecie 0 nie dostawałyby ruchów ani puszczenia przycisku
        on(events(fig).mousebutton; priority = 10_000) do ev
            ev.button == Makie.Mouse.left || return
            if ev.action == Makie.Mouse.press
                if in_axes()
                    x = mouse_X()
                    dragging[] = true
                    drag_x0[] = x
                    sel[] = (x, x)
                end
            elseif ev.action == Makie.Mouse.release
                dragging[] = false
                sel_vis[] && (sel_vis[] = false)
            end
        end
        on(events(fig).mouseposition; priority = 10_000) do _
            dragging[] || return
            if !(Makie.Mouse.left in events(fig).mousebuttonstate)    # puszczenie przegapione
                dragging[] = false
                sel_vis[] && (sel_vis[] = false)
                return
            end
            xa, xb = view_to_X(view_x[])
            sel[] = (drag_x0[], clamp(mouse_X(), xa, xb))
            sel_vis[] || (sel_vis[] = true)
        end

        # "Rozgrzej": wykonuje w podłączonej przeglądarce przykładowe operacje (przesuwanie, zoom, skala, jednostki,
        # jony, inny pierwiastek), żeby skompilować ścieżki kodu, a potem przywraca ustawienia. Działa raz na proces Julii.
        function warm_up!()
            warm_state[] == "idle" || return nothing
            warm_state[] = "running"
            saved = (view = view_x[], Z = Zsel[], xs = xs_d.idx[], unit = unit_d.idx[], ion2 = ion2_cb.value[],
                     T = T_sl.value[], thr = thr_sl.value[])
            step(f) = try
                f()
            catch e
                @warn "Rozgrzewanie: krok nie powiódł się" exception = (e, catch_backtrace())
            end
            try
                # przesunięcia ◀ ▶ oraz powiększenia i pomniejszenia ＋ －
                for (fac, sh) in ((1.0, 0.25), (1.0, 0.25), (1.0, -0.25), (0.5, 0.0), (0.5, 0.0), (2.0, 0.0), (2.0, 0.0), (1.0, -0.25))
                    step(() -> move_view!(fac, sh))
                end
                # zaznaczenie obszaru lewym przyciskiem = ustawienie limitów osi (obsługuje je ax.finallimits)
                for frac in (0.3, 0.3)
                    step() do
                        xa, xb = view_to_X(view_x[])
                        w = xb - xa
                        xlims!(ax, xa + frac * w, xb - frac * w)
                        request_redraw!()
                    end
                end
                # skala logarytmiczna, inne jednostki, drugi stopień jonizacji
                step(() -> (xs_d.idx[] = 2))
                step(() -> (xs_d.idx[] = saved.xs))
                step(() -> (unit_d.idx[] = 4))
                step(() -> (unit_d.idx[] = 5))
                step(() -> (unit_d.idx[] = saved.unit))
                step(() -> (ion2_cb.value[] = true))
                step(() -> (ion2_cb.value[] = saved.ion2))
                # suwak temperatury (przelicza natężenia LTE i udziały Sahy) oraz suwak zakresu jasności
                for T in (8000, 12000, 4000)
                    step(() -> (T_sl.value[] = T))
                end
                step(() -> (T_sl.value[] = saved.T))
                for thr in (3.0, 5.0, 2.0)
                    step(() -> (thr_sl.value[] = thr))
                end
                step(() -> (thr_sl.value[] = saved.thr))
                # bogaty pierwiastek (dużo linii: większa pojemność wykresu), potem powrót
                if 26 in ELEMENTS.Z
                    step(() -> set_element!(26))
                    step(() -> set_element!(saved.Z))
                end
            finally
                view_x[] = saved.view
                step(request_redraw!)
                WARMED[] = true
                warm_state[] = "done"
            end
            return nothing
        end
        for b in (warm_a, warm_b)
            on(b.value) do _
                @async warm_up!()
            end
        end

        # ---------------- tabela linii, kopiowanie, link ----------------
        csv_text(t; delim = ',') = (io = IOBuffer(); CSV.write(io, t; delim = delim); String(take!(io)))

        # Pobranie pliku CSV. Dane przekazujemy jako base64: Bonito wstawia interpolowane teksty do kodu JS w apostrofach
        # bez zamiany znaków nowej linii, więc wieloliniowy tekst bezpośrednio psułby składnię. Plik ma UTF-8 z BOM,
        # żeby Excel poprawnie czytał polskie i specjalne znaki.
        function download_text(text, fname; mime = "text/csv;charset=utf-8")
            b64 = Base64.base64encode(vcat(UInt8[0xef, 0xbb, 0xbf], codeunits(text)))
            Bonito.evaljs(session, Bonito.js"""
                (function () {
                    const bin = atob($(b64));
                    const bytes = new Uint8Array(bin.length);
                    for (let i = 0; i < bin.length; i++) { bytes[i] = bin.charCodeAt(i); }
                    const blob = new Blob([bytes], {type: $(mime)});
                    const url = URL.createObjectURL(blob);
                    const a = document.createElement('a');
                    a.href = url; a.download = $(fname);
                    document.body.appendChild(a); a.click(); a.remove();
                    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
                })();
            """)
        end

        for b in (table_a, table_b)
            on(b.value) do _
                table_open[] = !table_open[]
            end
        end
        for b in (tbl_close_a, tbl_close_b)
            on(b.value) do _
                table_open[] = false
            end
        end
        function export_table!(delim, ext, mime)       # pobranie tabeli linii z bieżącego widoku (wszystkie linie, nie tylko 200)
            t = current_table()
            t === nothing && return nothing
            x0, x1 = view_x[]
            sym = string(ELEMENTS.symbol[something(findfirst(==(Zsel[]), ELEMENTS.Z), 1)])
            download_text(csv_text(t; delim = delim), "spektrus_$(sym)_lines_$(round(Int, x0))-$(round(Int, x1))nm.$ext"; mime = mime)
            return nothing
        end
        for b in (tbl_csv_a, tbl_csv_b)
            on(b.value) do _
                export_table!(',', "csv", "text/csv;charset=utf-8")
            end
        end
        for b in (tbl_tsv_a, tbl_tsv_b)
            on(b.value) do _
                export_table!('\t', "tsv", "text/tab-separated-values;charset=utf-8")
            end
        end
        on(link_btn.value) do _
            qs = current_qs()
            url_last[] = qs
            Bonito.evaljs(session, Bonito.js"""
                (function () {
                    history.replaceState(null, '', '?' + $(qs));
                    const text = window.location.href;
                    function fallback() {
                        const ta = document.createElement('textarea');
                        ta.value = text; document.body.appendChild(ta); ta.select();
                        try { document.execCommand('copy'); } catch (e) {}
                        ta.remove();
                    }
                    if (navigator.clipboard && navigator.clipboard.writeText) {
                        navigator.clipboard.writeText(text).catch(fallback);
                    } else { fallback(); }
                })();
            """)
            link_msg[] = L("✓ skopiowano", "✓ copied")
            @async begin
                sleep(2.5)
                link_msg[] = ""
            end
        end

        # przykłady: ustawiają pierwiastek, stopnie jonizacji, światło widzialne i temperaturę, a resztę widoku na wartości domyślne
        function apply_preset!(Z, T, ions)
            preset_busy[] = true
            try
                set_element!(Z)
                T_sl.value[] = T
                ne_sl.value[] = 15.0
                saha_cb.value[] = true
                for (i, cb) in enumerate((ion1_cb, ion2_cb, ion3_cb, ion4_cb, ion5_cb))
                    cb.value[] = (i in ions)
                end
                unit_d.idx[] = 1;  med_d.idx[] = 1;  xs_d.idx[] = 1
                scale_d.idx[] = 1;  src_d.idx[] = 1;  ref_d.idx[] = 1
                thr_sl.value[] = 4.0
                view_x[] = RANGES[2][3]
            finally
                preset_busy[] = false
            end
            request_redraw!()
            return nothing
        end
        for (ba, bb, Z, T, ions) in preset_btns
            for b in (ba, bb)
                on(b.value) do _
                    Bonito.evaljs(session, Bonito.js"""document.querySelector('.sp-presets').removeAttribute('open');""")
                    @async apply_preset!(Z, T, ions)
                end
            end
        end

        # skróty klawiszowe: strzałki przesuwają, + i − zmieniają powiększenie, Home/End zakresy, F1 pomoc, Esc zamyka okna
        on(key_event) do s
            k = first(split(s, ':'))
            if k == "ArrowLeft"
                move_view!(1.0, -0.1)
            elseif k == "ArrowRight"
                move_view!(1.0, 0.1)
            elseif k == "+" || k == "="
                move_view!(0.5, 0.0)
            elseif k == "-" || k == "_"
                move_view!(2.0, 0.0)
            elseif k == "Home"
                view_x[] = RANGES[2][3]
                request_redraw!()
            elseif k == "End"
                view_x[] = RANGES[4][3]
                request_redraw!()
            elseif k == "F1"
                help_open[] = !help_open[]
                help_open[] && (about_open[] = false)
            elseif k == "Escape"
                about_open[] = false
                help_open[] = false
                table_open[] = false
            end
        end
        keys_js = Bonito.js"""
            (function () {
                document.addEventListener('keydown', function (e) {
                    const tag = ((e.target && e.target.tagName) || '').toLowerCase();
                    if (tag === 'input' || tag === 'textarea' || tag === 'select') { return; }
                    if (e.ctrlKey || e.metaKey || e.altKey) { return; }
                    const keys = ['ArrowLeft', 'ArrowRight', '+', '=', '-', '_', 'Home', 'End', 'F1', 'Escape'];
                    if (keys.indexOf(e.key) < 0) { return; }
                    e.preventDefault();
                    $(key_event).notify(e.key + ':' + Date.now());
                });
            })();
        """
        if isdefined(Bonito, :on_document_load)
            Bonito.on_document_load(session, keys_js)
        else
            @async (sleep(3); Bonito.evaljs(session, keys_js))
        end

        # zapis wykresu jako PNG: pobieramy obraz z canvasu WebGL w przeglądarce. WebGL czyści bufor rysowania po
        # wyświetleniu klatki, więc na chwilę podpinamy się pod requestAnimationFrame i robimy zrzut tuż po tym,
        # jak WGLMakie narysuje klatkę (pusty zrzut rozpoznajemy po rozmiarze i próbujemy w następnej klatce).
        # Wariant "jasny": na czas zrzutu wykres przełącza się na jasny motyw (białe tło, czarne osie, ciemniejsze
        # kreski), a po zrzucie (przeglądarka zwiększa licznik png_done) wraca do ciemnego.
        file_unit = Dict("nm" => "nm", "Å" => "A", "µm" => "um", "cm⁻¹" => "cm-1", "eV" => "eV")
        png_done  = Observable(0)
        export_token = Ref(0.0)

        function apply_theme!(light)
            fg = light ? :black : :white
            fig.scene.backgroundcolor[] = light ? RGBAf(1, 1, 1, 1) : RGBAf(0.09, 0.09, 0.11, 1)
            ax.backgroundcolor = light ? :white : :black
            ax_strip.backgroundcolor = light ? :white : :black
            ax_strip.leftspinecolor = light ? :gray30 : :gray40
            ax_strip.bottomspinecolor = light ? :gray30 : :gray40
            ax_strip.topspinecolor = light ? :gray30 : :gray40
            ax_strip.rightspinecolor = light ? :gray30 : :gray40
            ax.xlabelcolor = fg;      ax.ylabelcolor = fg
            ax.xticklabelcolor = fg;  ax.yticklabelcolor = fg
            ax.xtickcolor = fg;       ax.ytickcolor = fg
            ax.leftspinecolor = light ? :gray30 : :gray40
            ax.bottomspinecolor = light ? :gray30 : :gray40
            return nothing
        end
        function restore_theme!()
            print_mode[] || return nothing
            print_mode[] = false
            apply_theme!(false)
            request_redraw!()
            return nothing
        end
        on(png_done) do _
            restore_theme!()
        end

        function export_png!(light)
            Bonito.evaljs(session, Bonito.js"""document.querySelector('.sp-save').removeAttribute('open');""")   # zamknij menu
            sym = string(ELEMENTS.symbol[something(findfirst(==(Zsel[]), ELEMENTS.Z), 1)])
            fname = "spektrus_$(sym)_$(get(file_unit, unit_x[], "nm"))$(light ? "_light" : "").png"
            if light
                print_mode[] = true
                apply_theme!(true)
                request_redraw!()            # przerysuj kreski w kolorach do druku
                sleep(0.6)                   # daj przeglądarce czas na odebranie zmian
                token = time()
                export_token[] = token
                @async begin                 # zabezpieczenie: wróć do ciemnego motywu, gdyby zrzut się nie udał
                    sleep(10)
                    export_token[] == token && restore_theme!()
                end
            end
            Bonito.evaljs(session, Bonito.js"""
                (function () {
                    const canvas = document.querySelector('.sp-fig canvas');
                    const raf = window.requestAnimationFrame;
                    let frames = 0, done = false;
                    function restore() {
                        done = true;
                        window.requestAnimationFrame = raf;
                        $(png_done).notify(($(png_done).value || 0) + 1);      // zakończono (także po błędzie)
                    }
                    if (!canvas) { console.warn('PNG: nie znaleziono wykresu'); restore(); return; }
                    const blank = document.createElement('canvas');
                    blank.width = canvas.width; blank.height = canvas.height;
                    const blankLen = blank.toDataURL('image/png').length;
                    function save(url) {
                        restore();
                        const a = document.createElement('a');
                        a.href = url; a.download = $(fname);
                        document.body.appendChild(a); a.click(); a.remove();
                    }
                    window.requestAnimationFrame = function (cb) {
                        return raf.call(window, function (t) {
                            cb(t);
                            if (done) { return; }
                            frames += 1;
                            const url = canvas.toDataURL('image/png');
                            if (url.length > 1.3 * blankLen) { save(url); }
                            else if (frames > 90) { restore(); console.warn('PNG: pusty bufor WebGL, zrzut się nie udał'); }
                        });
                    };
                    setTimeout(function () { if (!done) { restore(); console.warn('PNG: przekroczono czas oczekiwania'); } }, 5000);
                })();
            """)
            return nothing
        end
        for b in (pngd_a, pngd_b)
            on(b.value) do _
                @async export_png!(false)
            end
        end
        for b in (pngl_a, pngl_b)
            on(b.value) do _
                @async export_png!(true)
            end
        end

        # zmiana języka: odśwież napisy generowane w Julii (opis linii, pasek informacji, etykiety osi)
        on(lang) do _
            set_hover!(L("Najedź kursorem na linię.", "Hover over a line."))
            set_summary!("", "")
            set_members!("")
            request_redraw!()
        end

        load_lap("definicje funkcji i reakcji")
        redraw!()    # pierwsze rysowanie
        load_lap("pierwsze rysowanie (redraw!)")
        url_ready[] = true

        # ---------------- układ strony ----------------
        # Cała strona mieści się w oknie (bez przewijania, bo kółko myszy służy do zoomu):
        # pasek narzędzi u góry, wykres zajmuje resztę wysokości, opis linii na dole.
        # Rzadko używane ustawienia są w rozwijanym panelu "⚙ Ustawienia" (nakłada się na wykres).
        css = Bonito.DOM.style("""
            html, body { margin: 0 !important; overflow: hidden !important; }
            .sp button { min-width: 0 !important; width: auto !important; padding: 3px 10px !important; margin: 0 !important; font-size: 0.9em !important; }
            .sp select { padding: 2px 4px !important; font-size: 0.9em !important; }
            .sp summary { cursor: pointer; user-select: none; padding: 3px 10px; border: 1px solid #666; border-radius: 4px; background: #2a2a31; }
            .sp summary:hover { background: #35353d; }
            .sp .sym select { width: 4.4em !important; min-width: 0 !important; }
            .sp .pair button, .sp .pair select { width: 100% !important; box-sizing: border-box !important; }
            .sp .view button { width: 2.4em !important; padding-left: 0 !important; padding-right: 0 !important; }
            .sp details[open] summary { background: #3a3a44; }
            """)
        row(items...; mb = 6) = Bonito.DOM.div(items...;
            style = "display:flex;flex-wrap:wrap;align-items:center;gap:6px 10px;margin-bottom:$(mb)px;")
        lab(key, items...) = Bonito.DOM.span(lbl(key), items...;                     # key: klucz w TR
            style = "display:inline-flex;align-items:center;gap:6px;")
        cblab(key, cb) = Bonito.DOM.label(cb, Bonito.DOM.span(tr(key));            # pole wyboru PRZED opisem, klikalny tekst
            style = "display:inline-flex;align-items:center;gap:6px;cursor:pointer;")
        ion(rom, cb, tipkey) = Bonito.DOM.span(Bonito.DOM.span(rom), cb; title = tr(tipkey),
            style = "display:inline-flex;align-items:center;gap:3px;")
        val(obs) = Bonito.DOM.span(map(string, obs); style = "min-width:3em;")
        grp(items...) = Bonito.DOM.span(items...; style = "display:inline-flex;align-items:center;gap:6px;flex-wrap:wrap;")
        muted(key) = Bonito.DOM.span(haskey(TR, key) ? lbl(key) : key;     # klucz w TR albo gotowy tekst (np. "│")
                                     style = "color:#9aa;font-size:0.85em;")

        more = Bonito.DOM.details(
            Bonito.DOM.summary(lbl("settings")),
            Bonito.DOM.div(
                row(Bonito.DOM.b(tr("g_axis")); mb = 4),
                row(lab("axscale", ddv(xs_d)),
                    Bonito.DOM.span(lbl("medium"), ddv(med_d);
                                    style = map(unit_d.idx) do i          # ukryty dla cm⁻¹ i eV (nie zależą od ośrodka)
                                        k = UNIT_KEYS[clamp(i, 1, length(UNIT_KEYS))]
                                        unit_is_energy(k) ? "display:none;" : "display:inline-flex;align-items:center;gap:6px;"
                                    end)),
                row(Bonito.DOM.b(tr("g_int")); mb = 4),
                row(lab("source", ddv(src_d)),
                    Bonito.DOM.span(cblab("saha", saha_cb), lab("lg_ne", ne_sl, val(ne_sl.value));
                                    style = "display:inline-flex;align-items:center;gap:14px;"),
                    lab("scale", ddv(scale_d)),
                    lab("reference", ddv(ref_d))),
                row(Bonito.DOM.b(tr("g_merge")); mb = 4),
                row(Bonito.DOM.span(lbl("threshold"), thr_sl, val(thr_sl.value); title = tr("thresholdtip"),
                                    style = "display:inline-flex;align-items:center;gap:6px;"),
                    lab("resol", res_sl, val(res_sl.value)),
                    Bonito.DOM.span(cblab("resolpx", rescb), px_sl, val(px_sl.value);
                                    style = "display:inline-flex;align-items:center;gap:6px;")),
                row(Bonito.DOM.b(tr("g_look")); mb = 4),
                row(cblab("bw", bw_cb),
                    cblab("gray", gray_cb),
                    cblab("members", members_cb); mb = 0);
                style = "position:absolute;left:12px;z-index:20;margin-top:4px;padding:10px 14px;" *
                        "background:#1e1e24;border:1px solid #666;border-radius:6px;box-sizing:border-box;" *
                        "width:max-content;max-width:calc(100vw - 24px);max-height:70vh;overflow:auto;" *
                        "box-shadow:0 6px 20px rgba(0,0,0,0.6);");
            style = "position:static;")

        help_panel = map(help_open, lang) do o, _l
            o || return Bonito.DOM.div()
            kbd(s) = Bonito.DOM.span(s; style = "border:1px solid #777;border-radius:3px;padding:0 5px;margin:0 1px;font-family:monospace;background:#2a2a31;")
            Bonito.DOM.div(
                Bonito.DOM.h3(L("Pomoc", "Help"); style = "margin:0 0 6px 0;"),
                Bonito.DOM.h4(L("Jak korzystać z programu", "How to use the program"); style = "margin:10px 0 4px 0;"),
                Bonito.DOM.ul(
                    Bonito.DOM.li(L("Wybierz pierwiastek. Nazwa, symbol i liczba atomowa są ze sobą powiązane, więc możesz użyć dowolnego z nich. " *
                                    "Zaznacz stopnie jonizacji: I to atomy obojętne, II–V to jony (II = atom, który stracił jeden elektron). " *
                                    "Zwykle wystarczy I.",
                                    "Choose an element. Name, symbol and atomic number are linked, so you can use any of them. " *
                                    "Tick the ionization stages: I is the neutral atom, II–V are ions (II = an atom that lost one electron). " *
                                    "Stage I is usually enough.")),
                    Bonito.DOM.li(L("Zakres widma: przyciski UV, widzialne, IR i wszystko. Widok przesuwasz przyciskami ◀ ▶, a powiększasz i pomniejszasz przyciskami ＋ i －. " *
                                    "Możesz też zaznaczyć fragment wykresu lewym przyciskiem myszy (zostanie powiększony), użyć kółka myszy do zoomu " *
                                    "i przeciągać prawym przyciskiem, żeby przesuwać.",
                                    "Spectrum range: the UV, visible, IR and all buttons. Pan with ◀ ▶ and zoom with ＋ and －. " *
                                    "You can also select a part of the plot with the left mouse button (it will be magnified), zoom with the mouse wheel " *
                                    "and drag with the right button to pan.")),
                    Bonito.DOM.li(L("Najedź kursorem na linię: pod wykresem pojawi się jej opis (pierwiastek i stopień, długość fali, przejście między poziomami, " *
                                    "współczynnik Aki i natężenie względne). Pod wykresem widać też położenie kursora w bieżącej jednostce.",
                                    "Hover over a line: its description appears below the plot (element and stage, wavelength, transition between levels, " *
                                    "Aki coefficient and relative intensity). The cursor position in the current unit is shown there as well.")),
                    Bonito.DOM.li(L("Linie leżące bardzo blisko siebie są rysowane jako jedna kreska. Powiększ widok albo zmniejsz rozdzielczość w ustawieniach, żeby je rozdzielić. " *
                                    "Dublet D sodu (589,0 i 589,6 nm) można tak rozdzielić.",
                                    "Lines lying very close together are drawn as one bar. Zoom in or lower the resolution in the settings to separate them. " *
                                    "The sodium D doublet (589.0 and 589.6 nm) can be separated this way.")),
                    Bonito.DOM.li(L("Jednostki osi: nm, Å, µm, cm⁻¹ i eV (dwie ostatnie rosną z energią, więc widmo jest odbite względem osi w nm).",
                                    "Axis units: nm, Å, µm, cm⁻¹ and eV (the last two increase with energy, so the spectrum is mirrored relative to the nm axis).")),
                    Bonito.DOM.li(L("⚙ Ustawienia: skala osi, powietrze lub próżnia, źródło natężeń, gęstość elektronowa n_e, zakres jasności " *
                                    "(ile rzędów wielkości słabsze od linii odniesienia linie pokazać), rozdzielczość zlewania linii i kolory. " *
                                    "Temperaturę T ustawiasz suwakiem obok.",
                                    "⚙ Settings: axis scale, air or vacuum, intensity source, electron density n_e, brightness range " *
                                    "(how many orders of magnitude weaker than the reference line to show), line merging resolution and colors. " *
                                    "The temperature T is set with the slider next to it.")),
                    Bonito.DOM.li(L("🧪 Przykłady jednym kliknięciem ustawiają lampę wodorową, helową, neonową, sodową albo rtęciową (stopień I) lub wapń z liniami H i K " *
                                    "Fraunhofera (jon Ca II, widoczny razem z atomem Ca I): pierwiastek, stopnie jonizacji, światło widzialne i poglądową " *
                                    "temperaturę. Przywracają też domyślne ustawienia widoku.",
                                    "🧪 Examples set up a hydrogen, helium, neon, sodium or mercury lamp (stage I) or calcium with the Fraunhofer H and K lines " *
                                    "(the Ca II ion, shown together with the Ca I atom) with one click: element, ionization stages, visible light and an " *
                                    "illustrative temperature. They also restore the default view settings.")),
                    Bonito.DOM.li(L("📋 Tabela linii pokazuje linie z widoku i pozwala pobrać je jako CSV lub TSV. 💾 Zapisz obraz zapisuje wykres jako PNG " *
                                    "(także na jasnym tle do druku). 🔗 Link kopiuje adres, który otwiera dokładnie ten sam widok.",
                                    "📋 Line table lists the lines in view and lets you download them as CSV or TSV. 💾 Save image saves the plot as a PNG " *
                                    "(also on a light background for printing). 🔗 Link copies an address that opens exactly the same view.")),
                    Bonito.DOM.li(L("Skróty klawiszowe: ", "Keyboard shortcuts: "),
                                  kbd("←"), kbd("→"), L(" przesuwanie (przytrzymaj, żeby przesuwać dalej); ", " pan (hold to keep panning); "),
                                  kbd("+"), kbd("−"), L(" powiększanie i pomniejszanie; ", " zoom in and out; "),
                                  kbd("Home"), L(" światło widzialne (380–780 nm); ", " visible light (380–780 nm); "),
                                  kbd("End"), L(" całe widmo; ", " whole spectrum; "),
                                  kbd("F1"), L(" ta pomoc; ", " this help; "),
                                  kbd("Esc"), L(" zamyka okna. Skróty nie działają, gdy kursor jest w polu tekstowym, na liście lub suwaku.",
                                                " closes the panels. Shortcuts do not work while a text field, list or slider has the focus.")),
                    style = "margin:0 0 6px 0;padding-left:20px;"),
                Bonito.DOM.h4(L("Krótkie wprowadzenie do spektroskopii", "A short introduction to spectroscopy"); style = "margin:12px 0 4px 0;"),
                Bonito.DOM.p(L("Światło białe rozszczepione przez pryzmat daje tęczę: widmo ciągłe, w którym kolory płynnie przechodzą jeden w drugi. " *
                               "Takie widmo wysyłają rozżarzone ciała stałe i ciecze oraz gęste gazy, na przykład włókno żarówki.",
                               "White light split by a prism gives a rainbow: a continuous spectrum in which the colors blend smoothly into one another. " *
                               "Glowing solids and liquids and dense gases, for example a light-bulb filament, emit such a spectrum.")),
                Bonito.DOM.p(L("Rozrzedzony gaz, który świeci (w wyładowaniu elektrycznym albo w płomieniu), wysyła inne światło: tylko kilka wąskich linii " *
                               "o różnych natężeniach, a pomiędzy nimi ciemność. To widmo liniowe, które pokazuje ten program. " *
                               "Przyrząd do oglądania widm to spektroskop; gdy zapisuje obraz, nazywa się spektrografem, a gdy mierzy natężenie " *
                               "poszczególnych części widma, spektrofotometrem.",
                               "A glowing low-density gas (in an electric discharge or a flame) emits different light: just a few narrow lines " *
                               "of different intensities with darkness between them. This is the line spectrum shown by this program. " *
                               "A device for viewing spectra is a spectroscope; if it records an image it is a spectrograph, and if it measures " *
                               "the intensity of different parts of the spectrum it is a spectrophotometer.")),
                Bonito.DOM.p(L("Linie powstają, bo energia atomu jest skwantowana: atom może mieć tylko ściśle określone energie. Gdy przechodzi z poziomu " *
                               "wyższego na niższy, wysyła foton o energii równej różnicy tych poziomów, czyli o ściśle określonej długości fali. " *
                               "Każdy pierwiastek ma inną budowę, więc inny zestaw linii: to jego odcisk palca. Dlatego z widma można rozpoznać skład " *
                               "chemiczny nieznanej substancji, także z daleka, na przykład gwiazdy (w widmach gwiazd linie najczęściej widać jako ciemne " *
                               "linie absorpcyjne w tych samych miejscach).",
                               "Lines arise because the energy of an atom is quantized: an atom can only have precisely defined energies. When it drops " *
                               "from a higher level to a lower one, it emits a photon whose energy equals the difference between the levels, so it has " *
                               "a precisely defined wavelength. Every element has a different structure and therefore a different set of lines: its " *
                               "fingerprint. That is why a spectrum can reveal the chemical composition of an unknown substance, even from afar, for " *
                               "example of a star (in stellar spectra the lines usually appear as dark absorption lines at the same wavelengths).")),
                Bonito.DOM.p(L("Spróbuj: wybierz sód. Najsilniejsza linia jest żółta (to w rzeczywistości dwie linie blisko siebie, powiększ widok), " *
                               "dlatego sodowe lampy uliczne świecą żółto-pomarańczowo. Rtęć ma silne linie w niebieskiej i zielonej części widma, " *
                               "a neon w pomarańczowo-czerwonej, stąd barwy lamp rtęciowych i neonowych.",
                               "Try it: choose sodium. The strongest line is yellow (really two lines close together, zoom in), which is why sodium " *
                               "street lamps glow yellow-orange. Mercury has strong lines in the blue and green parts of the spectrum and neon in the " *
                               "orange-red, hence the colors of mercury and neon lamps.")),
                Bonito.DOM.p(L("Światło nie kończy się na barwach, które widzimy. Program pokazuje też ultrafiolet i podczerwień; kolory tych linii są umowne " *
                               "(w ustawieniach można je zamienić na szare).",
                               "Light does not end at the colors we can see. The program also shows ultraviolet and infrared; the colors of those " *
                               "lines are conventional (you can turn them gray in the settings).")),
                Bonito.DOM.h4(L("O modelu natężeń", "About the intensity model"); style = "margin:12px 0 4px 0;"),
                Bonito.DOM.p(L("Natężenia linii liczy model równowagi termodynamicznej (LTE): zależą od temperatury T, od tego, ilu atomów jest w danym " *
                               "stopniu jonizacji (równowaga Saha, zależna od gęstości elektronowej n_e) i od prawdopodobieństwa przejścia Aki. " *
                               "To model poglądowy, dobry do nauki i porównań, a nie do pomiarów: w prawdziwych lampach i płomieniach natężenia mogą się " *
                               "różnić. Dla części pierwiastków i stopni NIST nie podaje Aki; wtedy używane są umowne natężenia z NIST i program " *
                               "ostrzega o tym pod wykresem.",
                               "Line intensities are computed with a local thermodynamic equilibrium (LTE) model: they depend on the temperature T, on how " *
                               "many atoms are in each ionization stage (Saha equilibrium, which depends on the electron density n_e) and on the " *
                               "transition probability Aki. It is an illustrative model, good for learning and comparisons but not for measurements: " *
                               "intensities in real lamps and flames can differ. For some elements and stages NIST gives no Aki; then the arbitrary " *
                               "NIST intensities are used and the program warns about it below the plot."));
                style = "padding:10px 14px;background:#1e1e24;border:1px solid #666;border-radius:6px;box-sizing:border-box;" *
                        "max-width:860px;max-height:70vh;overflow:auto;box-shadow:0 6px 20px rgba(0,0,0,0.6);line-height:1.4;")
        end

        table_panel = map(table_open, table_rev, lang) do o, _r, _l
            o || return Bonito.DOM.div()
            t = current_table()
            t === nothing && return Bonito.DOM.div(L("Brak linii w widoku.", "No lines in view."); style = "color:#aaa;")
            n = size(t, 1)
            m = min(n, TABLE_ROWS)
            unit = unit_x[]
            dig = get(Dict("nm" => 3, "Å" => 2, "µm" => 5, "cm⁻¹" => 2, "eV" => 4), unit, 3)
            th(s) = Bonito.DOM.th(s; style = "text-align:left;padding:2px 10px;border-bottom:1px solid #555;position:sticky;top:0;background:#1e1e24;")
            td(s, extra = "") = Bonito.DOM.td(s; style = "padding:1px 10px;white-space:nowrap;" * extra)
            aki(x) = ismissing(x) ? "—" : string(round(x; sigdigits = 3))
            rows = [Bonito.DOM.tr(td(string(round(t.axis_value[i]; digits = dig)), "font-family:monospace;text-align:right;"),
                                  td(fmt_rel(t.intensity_percent[i]), "text-align:right;"),
                                  td("$(t.element[i]) $(t.stage[i])"),
                                  td("$(t.lower_level[i]) → $(t.upper_level[i])", "color:#bbb;"),
                                  td(aki(t.Aki_s[i]), "text-align:right;color:#bbb;")) for i in 1:m]
            head = L("$n linii w widoku" * (n > m ? " (pokazano $m najsilniejszych; plik CSV zawiera wszystkie)" : ""),
                     "$n lines in view" * (n > m ? " ($m strongest shown; the CSV file contains all of them)" : ""))
            Bonito.DOM.div(
                Bonito.DOM.div(head; style = "color:#9aa;font-size:0.85em;margin-bottom:4px;"),
                Bonito.DOM.table(
                    Bonito.DOM.thead(Bonito.DOM.tr(th("λ [$unit]"), th("I [%]"), th(L("Linia", "Line")),
                                                   th(L("Przejście", "Transition")), th("Aki [s⁻¹]"))),
                    Bonito.DOM.tbody(rows...);
                    style = "border-collapse:collapse;font-size:0.85em;"))
        end
        table_box = Bonito.DOM.div(
            Bonito.DOM.div(Bonito.DOM.b(tr("tbl_title")),
                           pair(tbl_csv_a, tbl_csv_b; w = length("Download CSV") + 3),
                           pair(tbl_tsv_a, tbl_tsv_b; w = length("Download TSV") + 3),
                           pair(tbl_close_a, tbl_close_b; w = length("✕ Zamknij") + 3);
                           style = "display:flex;align-items:center;gap:10px;margin-bottom:6px;flex-wrap:wrap;"),
            Bonito.DOM.div(table_panel; style = "max-height:50vh;overflow:auto;");
            style = map(o -> o ? "display:block;background:#1e1e24;border:1px solid #666;border-radius:6px;padding:10px 14px;" *
                               "box-shadow:0 6px 20px rgba(0,0,0,0.6);box-sizing:border-box;" : "display:none;", table_open))

        w_preset = maximum(max(length(p[1]), length(p[2])) for p in PRESETS) + 3
        presets_menu = Bonito.DOM.details(
            Bonito.DOM.summary(lbl("presets")),
            Bonito.DOM.div(
                [row(pair(ba, bb; w = w_preset); mb = 4) for (ba, bb, _, _, _) in preset_btns]...,
                Bonito.DOM.div(tr("presetnote"); style = "color:#9aa;font-size:0.8em;max-width:380px;margin-top:6px;");
                style = "position:absolute;left:12px;z-index:20;margin-top:4px;padding:10px 14px;" *
                        "background:#1e1e24;border:1px solid #666;border-radius:6px;box-sizing:border-box;" *
                        "width:max-content;max-width:calc(100vw - 24px);box-shadow:0 6px 20px rgba(0,0,0,0.6);");
            class = "sp-presets", style = "position:static;")

        save_menu = Bonito.DOM.details(
            Bonito.DOM.summary(lbl("save")),
            Bonito.DOM.div(
                row(pair(pngd_a, pngd_b; w = length("Light background (for printing)") + 3); mb = 6),
                row(pair(pngl_a, pngl_b; w = length("Light background (for printing)") + 3); mb = 0);
                style = "position:absolute;left:12px;z-index:20;margin-top:4px;padding:10px 14px;" *
                        "background:#1e1e24;border:1px solid #666;border-radius:6px;box-sizing:border-box;" *
                        "width:max-content;max-width:calc(100vw - 24px);box-shadow:0 6px 20px rgba(0,0,0,0.6);");
            class = "sp-save", style = "position:static;")

        load_lap("panele i menu")
        return Bonito.DOM.div(
            css,
            row(Bonito.DOM.b("Spektruś.jl"; style = "font-size:1.2em;margin-right:6px;"),
                lab("element", ddv(name_d)),
                lab("symbol", Bonito.DOM.span(sym_dd; class = "sym")),
                Bonito.DOM.span(Bonito.DOM.span("Z:"), z_in; style = "display:inline-flex;align-items:center;gap:6px;"),
                grp(lbl("range"), range_nodes...),
                grp(lbl("view"), Bonito.DOM.span(b_left, b_right, b_zin, b_zout; class = "view", style = "display:inline-flex;gap:6px;")),
                lab("units", ddv(unit_d))),
            row(grp(lbl("lines"),
                    ion("I", ion1_cb, "tip1"),
                    muted("neutral"),
                    muted("│"),
                    muted("ionized"),
                    ion("II", ion2_cb, "tip2"),
                    ion("III", ion3_cb, "tip3"),
                    ion("IV", ion4_cb, "tip4"),
                    ion("V", ion5_cb, "tip5")),
                Bonito.DOM.span(Bonito.DOM.span("T [K]:"), T_sl, val(T_sl.value);
                                style = "display:inline-flex;align-items:center;gap:6px;"),
                more,
                presets_menu,
                pair(about_a, about_b; w = length("ⓘ O programie") + 3),
                pair(help_a, help_b; w = length("❓ Pomoc") + 3),
                Bonito.DOM.span(lang_btn; title = tr("langtip")),
                pair(table_a, table_b; w = length("📋 Line table") + 3),
                Bonito.DOM.span(link_btn; title = tr("linktip")),
                Bonito.DOM.span(link_msg; style = "color:#8f8;font-size:0.85em;"),
                save_menu,
                Bonito.DOM.span(pair(warm_a, warm_b; w = length("🔥 Rozgrzej") + 3); title = tr("warmtip"),
                                style = map(st -> st == "idle" ? "display:inline-block;" : "display:none;", warm_state)),
                Bonito.DOM.span(lbl("warming");
                                style = map(st -> st == "running" ? "display:inline-block;color:#fc8;" : "display:none;", warm_state))),
            Bonito.DOM.div(Bonito.DOM.div(about_panel; style = "position:absolute;z-index:30;left:12px;right:12px;");
                           style = "position:relative;height:0;"),
            Bonito.DOM.div(Bonito.DOM.div(help_panel; style = "position:absolute;z-index:31;left:12px;right:12px;");
                           style = "position:relative;height:0;"),
            Bonito.DOM.div(Bonito.DOM.div(table_box; style = "position:absolute;z-index:25;left:12px;right:12px;");
                           style = "position:relative;height:0;"),
            Bonito.DOM.div(fig; class = "sp-fig", style = "flex:1 1 auto;min-height:0;position:relative;overflow:hidden;"),
            Bonito.DOM.div(info; style = "margin-top:4px;color:#9aa;font-size:0.85em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;"),
            Bonito.DOM.div(note; style = "color:#fc8;font-size:0.85em;overflow:hidden;"),
            Bonito.DOM.div(
                Bonito.DOM.div(Bonito.DOM.span(cursor; style = "float:right;color:#9aa;font-size:0.85em;font-family:monospace;"),
                               summary; style = "min-height:1.3em;"),
                Bonito.DOM.div(hover; style = "font-family:monospace;color:#bbb;min-height:1.3em;"),
                Bonito.DOM.div(members_txt; style = "font-family:monospace;color:#bbb;white-space:pre;overflow:auto;max-height:6em;"),
                Bonito.DOM.div(hint; style = "color:#889;font-size:0.9em;min-height:1.2em;");
                style = "margin-top:4px;padding:4px 10px;border:1px solid #444;border-radius:6px;font-size:0.92em;flex:none;max-height:30vh;overflow:auto;");
            class = "sp",
            style = "position:relative;display:flex;flex-direction:column;height:100vh;box-sizing:border-box;overflow:hidden;" *
                    "background:#161618;color:#ddd;font-family:system-ui,sans-serif;padding:8px 12px;")
    end
end

# Wstępne przygotowanie rdzenia: wczytanie levels.csv i wykonanie głównych obliczeń na próbce, żeby kompilacja i
# wczytywanie danych nie obciążały pierwszej strony po starcie (część JIT dotyczy tylko Makie/Bonito i wymaga
# prawdziwej strony, tego nie da się rozgrzać po stronie serwera).
function warm_core!()
    t0 = time()
    try
        level_db()
        for (Z, ions) in ((11, [1]), (26, [1, 2]))
            sp = prepare_lines(DF, Z; ions = ions, medium = "air", source = "lte", T = 6000.0, saha = true, ne = 1e15)
            l = sp.lines
            isempty(l.I) && continue
            cl = cluster_window(l, -Inf, Inf, 0.1)
            Iref = reference_intensity(l, 380.0, 780.0; ref = "visible")
            idx = findall(visible_mask(cl.I, Iref, 4.0))
            isempty(idx) || line_height(cl.I[first(idx)], Iref, 4.0; scale = "log")
            line_table(l, collect(1:min(100, size(l, 1))), Iref, "nm")
        end
        fig = Figure(size = (600, 300))                  # kompilacja tworzenia wykresów Makie (bez przeglądarki)
        ax = Axis(fig[1, 1])
        linesegments!(ax, [Point2d(0, 0), Point2d(1, 1)]; color = [RGBAf(1, 1, 1, 1), RGBAf(1, 1, 1, 1)])
        xlims!(ax, 0, 1)
    catch e
        @warn "Wstępne przygotowanie nie powiodło się (aplikacja zadziała, tylko wolniej na starcie)" exception = (e, catch_backtrace())
    end
    println("Przygotowanie rdzenia: $(round(time() - t0; digits = 1)) s")
end

app = build_app()

if abspath(PROGRAM_FILE) == @__FILE__
    get(ENV, "SPEKTRUS_WARM", "1") == "0" || warm_core!()        # SPEKTRUS_WARM=0 wyłącza przygotowanie rdzenia przy starcie
    server = Bonito.Server(app, "127.0.0.1", 8080)
    println("Spektruś działa: http://127.0.0.1:8080   (Ctrl+C kończy)")
    try
        while true
            sleep(1)
        end
    catch e
        e isa InterruptException || rethrow()
    end
end