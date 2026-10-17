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
const RES_FLOOR   = 1e-4      # najmniejsza rozdzielczość zlewania linii [nm] (linie bliższe niż to uznajemy za tę samą)

# autor tej wersji (wyświetlany w oknie "O programie" obok Claude'a) - wpisz swoje imię i nazwisko
const AUTHOR = "WPISZ IMIĘ I NAZWISKO"

const DF       = load_data()
const ELEMENTS = element_list(DF)
const NAME_OPTS = [element_name(Z) for Z in ELEMENTS.Z]
const SYM_OPTS  = collect(String, ELEMENTS.symbol)
const RANGES   = ["UV" => (100.0, 380.0), "widzialne" => (380.0, 780.0),
                  "IR" => (780.0, 2500.0), "wszystko" => (100.0, 2500.0)]

# jednostki osi X: klucze (używane w kodzie) i napisy na liście
const UNIT_KEYS = ["nm", "Å", "µm", "cm⁻¹", "eV"]
const UNIT_OPTS = ["nanometry (nm)", "angstremy (Å)", "mikrometry (µm)", "liczba falowa (cm⁻¹)", "energia fotonu (eV)"]

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

# „1 linia spektralna”, „2 linie spektralne”, „5 linii spektralnych”
function lines_pl(n)
    n == 1 && return "1 linia spektralna"
    (2 <= n % 10 <= 4 && !(12 <= n % 100 <= 14)) && return "$n linie spektralne"
    return "$n linii spektralnych"
end

# względna suma natężeń: liczba całkowita od 10 wzwyż, poniżej 2 cyfry znaczące
fmt_rel(v) = !isfinite(v) ? "∞" : v >= 1e6 ? string(round(v; sigdigits = 3)) :
             v >= 10 ? string(round(Int, v)) : string(round(v; sigdigits = 2))

const ROMAN = ("I", "II", "III", "IV", "V")

# udział procentowy: pełne procenty od 10% wzwyż, poniżej 2 cyfry znaczące
fmt_pct(f) = f >= 0.1 ? "$(round(Int, 100f))%" : "$(round(100f; sigdigits = 2))%"

function build_app()
    return Bonito.App() do session
        # ---------------- kontrolki ----------------
        # pierwiastek: nazwa, symbol i liczba atomowa to trzy powiązane kontrolki (jedna wspólna wartość Zsel)
        z0       = something(findfirst(==(11), ELEMENTS.Z), 1)
        Zsel     = Observable(ELEMENTS.Z[z0])
        name_dd  = Bonito.Dropdown(NAME_OPTS; index = z0)
        sym_dd   = Bonito.Dropdown(SYM_OPTS; index = z0)
        # liczba atomowa: własne pole <input type=number> z Observable{Int} (Bonito.NumberInput istnieje tylko dla Float64)
        z_val    = Observable(Int(ELEMENTS.Z[z0]))
        z_in     = Bonito.DOM.input(; type = "number", value = z_val, min = 1, max = 99, step = 1,
                                    style = "width:3.5em;padding:2px 4px;font-weight:600;font-size:0.9em;border-radius:3px;",
                                    onchange = Bonito.js"""event => {
                                        const v = Math.round(parseFloat(event.srcElement.value));
                                        if (Number.isFinite(v) && $(z_val).value != v) { $(z_val).notify(v); }
                                    }""")
        unit_dd  = Bonito.Dropdown(UNIT_OPTS)
        med_dd   = Bonito.Dropdown(["powietrze", "próżnia"])
        xs_dd    = Bonito.Dropdown(["liniowa", "logarytmiczna"])
        scale_dd = Bonito.Dropdown(["logarytmiczna", "liniowa"])
        src_dd   = Bonito.Dropdown(["model LTE", "NIST"])
        ref_dd   = Bonito.Dropdown(["linia widzialna", "linia w oknie", "linia pierwiastka"])
        ion1_cb  = Bonito.Checkbox(true)                     # linie niezjonizowane (I)
        ion2_cb  = Bonito.Checkbox(false)                    # linie zjonizowane (II)
        ion3_cb  = Bonito.Checkbox(false)                    # podwójnie zjonizowane (III)
        ion4_cb  = Bonito.Checkbox(false)                    # potrójnie zjonizowane (IV)
        ion5_cb  = Bonito.Checkbox(false)                    # poczwórnie zjonizowane (V)
        saha_cb  = Bonito.Checkbox(true)                     # równowaga jonizacji (Saha)
        bw_cb    = Bonito.Checkbox(false)                    # tryb czarno-biały
        members_cb = Bonito.Checkbox(false)                  # lista składowych zlanych linii w opisie
        T_sl     = Bonito.Slider(1000:100:20000; value = 6000)
        ne_sl    = Bonito.Slider(12.0:0.5:18.0; value = 15.0)     # lg n_e [cm⁻³] (dla Sahy)
        res_sl   = Bonito.Slider(0.0:0.1:5.0; value = 0.0)
        thr_sl   = Bonito.Slider(0.5:0.1:6.0; value = 4.0)
        rescb    = Bonito.Checkbox(false)                    # rozdzielczość zależna od zoomu
        px_sl    = Bonito.Slider(1:10; value = 2)            # ... w pikselach ekranu
        btns     = [Bonito.Button(name) for (name, _) in RANGES]
        b_left   = Bonito.Button("◀")
        b_right  = Bonito.Button("▶")
        b_zin    = Bonito.Button("＋")
        b_zout   = Bonito.Button("－")
        about_btn = Bonito.Button("ⓘ O programie")
        about_open = Observable(false)

        info  = Observable("")
        hover = Observable("Najedź kursorem na linię.")
        summary = Observable("")      # „N linii spektralnych o względnej sumie natężeń X”
        hint    = Observable("")      # podpowiedź o zaznaczaniu i powiększaniu
        members_txt = Observable("")  # składowe zlanych linii (gdy włączone)

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

        # ---------------- stan ----------------
        view_x    = Ref((380.0, 780.0))   # bieżący zakres osi X [nm] (w nm, nie w log10)
        redrawing = Ref(false)            # ignorujemy zmiany limitów wywołane przez nasz kod
        cache     = Ref{Any}(nothing)     # ostatnio policzone linie (element, T, ośrodek, źródło)
        shown     = Ref{Any}(nothing)     # to, co widać na wykresie (dla podpowiedzi)
        cur_res   = Ref(RES_FLOOR)        # rozdzielczość zlewania użyta przy ostatnim rysowaniu [nm]

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
            spread == 0 && return "Te linie mają identyczną długość fali w danych, więc nie da się ich rozdzielić."
            cur_res[] <= RES_FLOOR * 1.001 && return ""          # najdrobniejsza rozdzielczość: nic już nie rozdzieli
            rescb.value[] && return "Możesz zaznaczyć i powiększyć ten wycinek w celu zwiększenia rozdzielczości."
            return "Zmniejsz rozdzielczość (suwak) albo włącz jej zależność od zoomu, żeby rozdzielić te linie."
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
                medium = med_dd.value[] == "powietrze" ? "air" : "vac"
                source = src_dd.value[] == "model LTE" ? "lte" : "nist"
                scale  = scale_dd.value[] == "logarytmiczna" ? "log" : "lin"
                refmode = ref_dd.value[] == "linia w oknie" ? "window" :
                          ref_dd.value[] == "linia pierwiastka" ? "all" : "visible"
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
                logmode = xs_dd.value[] == "logarytmiczna"
                unit_key = UNIT_KEYS[something(findfirst(==(unit_dd.value[]), UNIT_OPTS), 1)]
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
                    (cr, cg, cb) = bw ? (1.0, 1.0, 1.0) : wl2rgb(cl.λ[k])
                    a = 0.15 + 0.85 * h
                    i2 = 2 * j
                    pts1[i2 - 1] = Point2d(X, 0.0);  pts1[i2] = Point2d(X, 1.0)
                    col1[i2 - 1] = RGBAf(cr, cg, cb, a);  col1[i2] = RGBAf(cr, cg, cb, a)
                    pts2[i2 - 1] = Point2d(X, 0.0);  pts2[i2] = Point2d(X, h)
                    col2[i2 - 1] = RGBAf(cr, cg, cb, 1);  col2[i2] = RGBAf(cr, cg, cb, 1)
                end
                lap("punkty i kolory ($nseg kresek)")
                sg.pts1[] = pts1;  sg.col1[] = col1
                sg.pts2[] = pts2;  sg.col2[] = col2
                lap("aktualizacja danych wykresów")
                xlabel_txt = unit_key == "cm⁻¹" ? "liczba falowa σ = 1/λ [cm⁻¹]" :
                             unit_key == "eV" ? "energia fotonu [eV]" :
                             "λ [$unit_key, $(medium == "air" ? "powietrze" : "próżnia")]"
                ax.xlabel[] == xlabel_txt || (ax.xlabel = xlabel_txt)
                xa, xb = view_to_X((x0, x1))
                xlims!(ax, xa, xb)
                ax.xticks = ticks_for(x0, x1)
                ylims!(ax, 0, 1.05)
                ylims!(ax_strip, 0, 1)
                lap("limity i ticki")

                shown[] = (λ = cl.λ[idx], I = cl.I[idx], Iref = Iref, rows = cl.rows[idx], n = cl.n[idx], spread = cl.spread[idx],
                           first = cl.first[idx], last = cl.last[idx], lines = lines)
                frac = cache[].spec.fractions
                saha_txt = frac === nothing ? "" :
                    ", Saha: lg n_e = $(round(ne_exp; digits = 1)), udziały " *
                    join(["$(ROMAN[z]) $(fmt_pct(frac[z]))" for z in eachindex(frac) if frac[z] >= 0.0005], ", ")
                src_txt = isempty(lines.I) ? "brak linii (dla wybranych stopni jonizacji nie ma danych)" :
                          use_lte ? "model LTE, T = $(Int(T)) K$saha_txt" :
                          (source == "lte" ? "NIST (brak Aki dla tego pierwiastka, LTE niedostępny)" : "NIST")
                info[] = "Linie widoczne: $(length(idx)) z $(length(cl.λ)) (całe widmo)   ·   " *
                         "źródło natężeń: $src_txt   ·   odniesienie jasności: $(round(Iref; sigdigits = 3))   ·   " *
                         "rozdzielczość: $(round(res; sigdigits = 2)) nm"
                lap("info")
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
        onany(Zsel, unit_dd.value, med_dd.value, xs_dd.value, scale_dd.value, src_dd.value, ref_dd.value,
              ion1_cb.value, ion2_cb.value, ion3_cb.value, ion4_cb.value, ion5_cb.value,
              saha_cb.value, ne_sl.value, bw_cb.value, T_sl.value, res_sl.value, thr_sl.value, rescb.value, px_sl.value) do args...
            request_redraw!()
        end

        # powiązane kontrolki pierwiastka: nazwa, symbol i liczba atomowa zmieniają się razem
        # (jedna wspólna wartość Zsel; sprawdzanie równości zapobiega pętlom)
        function set_element!(Z)
            Z == Zsel[] || (Zsel[] = Z)
            return nothing
        end
        on(name_dd.option_index) do i
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
            name_dd.option_index[] == i || (name_dd.option_index[] = i)
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
        on(about_btn.value) do _
            about_open[] = !about_open[]
        end
        about_panel = map(about_open) do o
            o || return Bonito.DOM.div()
            Bonito.DOM.div(
                Bonito.DOM.h3("Spektruś.jl"; style = "margin:0 0 6px 0;"),
                Bonito.DOM.p(Bonito.DOM.b("Autorzy tej wersji: "), AUTHOR, " oraz Claude (Anthropic)."),
                Bonito.DOM.p("Odtworzenie programu Spektruś 1.0 (Mikołaj Pytel, 2002) w Julii. " *
                             "Pokazuje widma emisyjne pierwiastków 1–99 w UV, świetle widzialnym i podczerwieni."),
                Bonito.DOM.p("Dane linii spektralnych: ", Bonito.DOM.b("NIST Atomic Spectra Database"),
                             " (stopnie jonizacji I–V, 100–2500 nm). Prosimy cytować: A. Kramida, Yu. Ralchenko, J. Reader " *
                             "i NIST ASD Team, NIST Atomic Spectra Database, ",
                             Bonito.DOM.a("https://physics.nist.gov/asd"; href = "https://physics.nist.gov/asd",
                                          target = "_blank", style = "color:#8cf;"), "."),
                Bonito.DOM.p("Natężenia: model LTE, emisja optycznie cienka, I ∝ (N_z/N)·g_k·A_ki·ν·exp(−E_k/kT)/U_z(T) " *
                             "(temperatura z suwaka T). Udziały stopni jonizacji N_z/N wynikają z równowagi Saha dla " *
                             "gęstości elektronowej n_e (suwak lg n_e); uwzględniono stopnie I–V. Funkcje podziału U_z liczone " *
                             "są z poziomów występujących w pobranych liniach, więc to przybliżenie (dla złożonych atomów " *
                             "może zaniżać U). Energie jonizacji: tabela z pakietu mendeleev (zgodna z NIST ASD w sprawdzonych " *
                             "przypadkach). Po wyłączeniu opcji Saha jasności linii różnych stopni jonizacji nie są " *
                             "porównywalne, a w obrębie jednego stopnia są."),
                Bonito.DOM.p(Bonito.DOM.b("Źródła modelu fizycznego. "),
                             "Rozkład Boltzmanna populacji poziomów, emisja optycznie cienka w LTE oraz równowaga jonizacji Saha " *
                             "to standardowe elementy fizyki plazmy; opisano je m.in. w:"),
                Bonito.DOM.ul(
                    Bonito.DOM.li("H. R. Griem, ", Bonito.DOM.i("Principles of Plasma Spectroscopy"),
                                  ", Cambridge University Press, 1997 (rozdz. 7: Thermodynamic equilibrium relations), ",
                                  Bonito.DOM.a("strona książki";
                                               href = "https://www.cambridge.org/core/books/principles-of-plasma-spectroscopy/90F6DBB512089B7E3AC121744EFE3D93",
                                               target = "_blank", style = "color:#8cf;"), "."),
                    Bonito.DOM.li("M. N. Saha, ", Bonito.DOM.i("Ionization in the solar chromosphere"),
                                  ", Philosophical Magazine 40, 472 (1920), doi: ",
                                  Bonito.DOM.a("10.1080/14786441008636148";
                                               href = "https://doi.org/10.1080/14786441008636148",
                                               target = "_blank", style = "color:#8cf;"), ".");
                    style = "margin:4px 0 8px 0;padding-left:20px;"),
                Bonito.DOM.p("Stała w równaniu Saha: 2(2π m_e k / h²)^{3/2} = 4,8294·10¹⁵ cm⁻³ K⁻³ᐟ² (stałe CODATA). " *
                             "Model jest orientacyjny: zakłada LTE i optycznie cienką plazmę, które nie zawsze są spełnione."),
                Bonito.DOM.p("Obsługa: kółko myszy = zoom, przeciągnięcie lewym = zaznaczenie obszaru, " *
                             "prawym = przesuwanie; przyciski ◀ ▶ ＋ － działają bez myszy. Jednostki osi: nm, Å, µm, cm⁻¹, eV " *
                             "(cm⁻¹ i eV nie zależą od ośrodka: liczone są z długości fali w próżni, a oś rośnie z energią; " *
                             "przełącznik powietrze/próżnia jest wtedy ukryty).");
                style = "max-width:760px;margin:8px 0;padding:10px 14px;background:#222228;" *
                        "border:1px solid #444;border-radius:6px;font-size:0.92em;line-height:1.4;" *
                        "box-sizing:border-box;max-height:70vh;overflow:auto;")
        end

        for (b, (_, (lo, hi))) in zip(btns, RANGES)
            on(b.value) do _
                view_x[] = (lo, hi)
                request_redraw!()
            end
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
            ax.xticks = ticks_for(clo, chi)
            (ref_dd.value[] == "linia w oknie" || rescb.value[]) && schedule_redraw!()
        end

        # podpowiedź po najechaniu na linię
        function hover_handler(_)
            s = shown[]
            s === nothing && return
            if !(Makie.is_mouseinside(ax.scene) || Makie.is_mouseinside(ax_strip.scene))
                return
            end
            λs = s.λ
            if isempty(λs)
                set_hover!("Brak widocznych linii.")
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
                set_hover!(describe_line(s.lines, s.rows[i]; n = s.n[i], spread = s.spread[i], unit = unit_x[]))
                set_summary!("$(lines_pl(s.n[i])) o względnej sumie natężeń $(fmt_rel(100 * s.I[i] / s.Iref))",
                             merge_hint(s.n[i], s.spread[i]))
                set_members!(members_cb.value[] && s.n[i] > 1 ? describe_members(s.lines, s.first[i], s.last[i]) : "")
            else
                set_hover!("Najedź kursorem na linię.")
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

        redraw!()    # pierwsze rysowanie

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
            .sp details[open] summary { background: #3a3a44; }
            """)
        row(items...; mb = 6) = Bonito.DOM.div(items...;
            style = "display:flex;flex-wrap:wrap;align-items:center;gap:6px 14px;margin-bottom:$(mb)px;")
        lab(text, items...) = Bonito.DOM.span(Bonito.DOM.span(text * " "), items...;
            style = "display:inline-flex;align-items:center;gap:6px;")
        ion(rom, cb, tip) = Bonito.DOM.span(Bonito.DOM.span(rom), cb; title = tip,
            style = "display:inline-flex;align-items:center;gap:3px;")
        val(obs) = Bonito.DOM.span(map(string, obs); style = "min-width:3em;")
        grp(items...) = Bonito.DOM.span(items...; style = "display:inline-flex;align-items:center;gap:6px;flex-wrap:wrap;")
        muted(text) = Bonito.DOM.span(text; style = "color:#9aa;font-size:0.85em;")

        more = Bonito.DOM.details(
            Bonito.DOM.summary("⚙ Ustawienia"),
            Bonito.DOM.div(
                row(Bonito.DOM.b("Oś"); mb = 4),
                row(lab("Skala osi:", xs_dd),
                    Bonito.DOM.span(Bonito.DOM.span("Długość fali: "), med_dd;
                                    style = map(unit_dd.value) do v          # ukryty dla cm⁻¹ i eV (nie zależą od ośrodka)
                                        k = UNIT_KEYS[something(findfirst(==(v), UNIT_OPTS), 1)]
                                        unit_is_energy(k) ? "display:none;" : "display:inline-flex;align-items:center;gap:6px;"
                                    end)),
                row(Bonito.DOM.b("Natężenia"); mb = 4),
                row(lab("Źródło:", src_dd),
                    lab("Saha, lg n_e [cm⁻³]:", saha_cb, ne_sl, val(ne_sl.value)),
                    lab("Skala:", scale_dd),
                    lab("Odniesienie:", ref_dd)),
                row(Bonito.DOM.b("Zlewanie linii"); mb = 4),
                row(lab("Próg [dekady]:", thr_sl, val(thr_sl.value)),
                    lab("Rozdzielczość [nm]:", res_sl, val(res_sl.value)),
                    lab("lub zależna od zoomu [px]:", rescb, px_sl, val(px_sl.value))),
                row(Bonito.DOM.b("Wygląd"); mb = 4),
                row(lab("czarno-białe", bw_cb),
                    lab("składowe zlanych linii w opisie", members_cb); mb = 0);
                style = "position:absolute;left:12px;z-index:20;margin-top:4px;padding:10px 14px;" *
                        "background:#1e1e24;border:1px solid #666;border-radius:6px;box-sizing:border-box;" *
                        "width:max-content;max-width:calc(100vw - 24px);max-height:70vh;overflow:auto;" *
                        "box-shadow:0 6px 20px rgba(0,0,0,0.6);");
            style = "position:static;")

        return Bonito.DOM.div(
            css,
            row(Bonito.DOM.b("Spektruś.jl"; style = "font-size:1.2em;margin-right:6px;"),
                lab("Pierwiastek:", name_dd),
                lab("Symbol:", sym_dd),
                lab("Z:", z_in),
                grp(Bonito.DOM.span("Zakres:"), btns...),
                grp(Bonito.DOM.span("Widok:"), b_left, b_right, b_zin, b_zout),
                lab("Jednostki:", unit_dd)),
            row(grp(Bonito.DOM.span("Linie:"),
                    ion("I", ion1_cb, "atom niezjonizowany (X)"),
                    muted("niezjonizowane"),
                    muted("│"),
                    muted("zjonizowane:"),
                    ion("II", ion2_cb, "jednokrotnie zjonizowane (X⁺)"),
                    ion("III", ion3_cb, "dwukrotnie zjonizowane (X²⁺)"),
                    ion("IV", ion4_cb, "trzykrotnie zjonizowane (X³⁺)"),
                    ion("V", ion5_cb, "czterokrotnie zjonizowane (X⁴⁺)")),
                lab("T [K]:", T_sl, val(T_sl.value)),
                more,
                about_btn),
            Bonito.DOM.div(Bonito.DOM.div(about_panel; style = "position:absolute;z-index:30;left:12px;right:12px;");
                           style = "position:relative;height:0;"),
            Bonito.DOM.div(fig; style = "flex:1 1 auto;min-height:0;position:relative;overflow:hidden;"),
            Bonito.DOM.div(info; style = "margin-top:4px;color:#9aa;font-size:0.85em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;"),
            Bonito.DOM.div(
                Bonito.DOM.div(summary; style = "min-height:1.3em;"),
                Bonito.DOM.div(hover; style = "font-family:monospace;color:#bbb;min-height:1.3em;"),
                Bonito.DOM.div(members_txt; style = "font-family:monospace;color:#bbb;white-space:pre;overflow:auto;max-height:6em;"),
                Bonito.DOM.div(hint; style = "color:#889;font-size:0.9em;min-height:1.2em;");
                style = "margin-top:4px;padding:4px 10px;border:1px solid #444;border-radius:6px;font-size:0.92em;flex:none;max-height:30vh;overflow:auto;");
            class = "sp",
            style = "position:relative;display:flex;flex-direction:column;height:100vh;box-sizing:border-box;overflow:hidden;" *
                    "background:#161618;color:#ddd;font-family:system-ui,sans-serif;padding:8px 12px;")
    end
end

app = build_app()

if abspath(PROGRAM_FILE) == @__FILE__
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