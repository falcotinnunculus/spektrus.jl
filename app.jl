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

WGLMakie.activate!()

# ---- granice osi λ (zmień tutaj) ----
const X_MIN       = 0.0       # dolna granica [nm], skala liniowa
const X_MIN_LOG   = 10.0      # dolna granica [nm], skala logarytmiczna (log(0) nie istnieje)
const X_MAX       = 3000.0    # górna granica [nm]; dane NIST sięgają 2500 nm
const X_MIN_WIDTH = 0.002     # najwęższe okno [nm] w skali liniowej (w logarytmicznej: 0,002 dekady)
const RES_FLOOR   = 1e-4      # najmniejsza rozdzielczość zlewania linii [nm] (linie bliższe niż to uznajemy za tę samą)

const DF       = load_data()
const ELEMENTS = element_list(DF)
const OPTS     = ["$(r.symbol) ($(r.Z))" for r in eachrow(ELEMENTS)]
const RANGES   = ["UV" => (100.0, 380.0), "widzialne" => (380.0, 780.0),
                  "IR" => (780.0, 2500.0), "wszystko" => (100.0, 2500.0)]

"""
Przycina okno [lo, hi] [nm] do granic osi, zachowując jego szerokość (przy przesuwaniu zatrzymuje się
na brzegu). W skali logarytmicznej liczy w log10, więc zachowuje stosunek hi/lo.
"""
function clamp_view(lo, hi, logmode)
    t  = logmode ? log10 : identity
    ti = logmode ? (x -> 10.0^x) : identity
    if logmode
        lo = max(lo, X_MIN_LOG)
        hi = max(hi, lo * 1.001)
    end
    amin = t(logmode ? X_MIN_LOG : X_MIN)
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
        elem_dd  = Bonito.Dropdown(OPTS; index = something(findfirst(==(11), ELEMENTS.Z), 1))
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
                  xtickformat = fmt_nm, axstyle...)
        linkxaxes!(ax_strip, ax)
        xlims!(ax, 380, 780)

        # ---------------- stan ----------------
        view_x    = Ref((380.0, 780.0))   # bieżący zakres osi X [nm] (w nm, nie w log10)
        log_x     = Ref(false)            # czy oś X jest logarytmiczna
        redrawing = Ref(false)            # ignorujemy zmiany limitów wywołane przez nasz kod
        cache     = Ref{Any}(nothing)     # ostatnio policzone linie (element, T, ośrodek, źródło)
        shown     = Ref{Any}(nothing)     # to, co widać na wykresie (dla podpowiedzi)

        # Makie trzyma limity w przestrzeni przetransformowanej (log10 w skali logarytmicznej);
        # zabezpieczenie: wartości < 4 w skali log traktujemy jako log10(nm), bo dolna granica to 10 nm.
        function limits_to_nm(r)
            lo, hi = r.origin[1], r.origin[1] + r.widths[1]
            if log_x[] && hi < 4.0
                return (10.0^lo, 10.0^hi)
            end
            return (lo, hi)
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
        # zerowej długości w x = 1 nm (dodatnie, więc bezpieczne w skali log). Dzięki temu tablice punktów
        # i kolorów mają zawsze tę samą długość. Wykresy budujemy od nowa tylko przy zmianie skali osi X
        # albo gdy zabraknie pojemności.
        segs  = Ref{Any}(nothing)
        PAD   = Point2f(1, 0)
        BLANK = RGBAf(0, 0, 0, 0)

        # nakładka pokazująca zaznaczany obszar (przeciąganie lewym przyciskiem); x w nm,
        # wartości początkowe dodatnie, żeby nie psuły skali log
        sel       = Observable((1.0, 1.0))
        sel_vis   = Observable(false)
        sel_main  = lift(sel) do r
            Rect2f(Point2f(min(r[1], r[2]), 0), Vec2f(abs(r[2] - r[1]), 1.05))
        end
        sel_strip = lift(sel) do r
            Rect2f(Point2f(min(r[1], r[2]), 0), Vec2f(abs(r[2] - r[1]), 1))
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
        TIMING = true
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
        function set_summary!(sum_txt, hint_txt)
            summary[] == sum_txt || (summary[] = sum_txt)
            hint[] == hint_txt || (hint[] = hint_txt)
            return nothing
        end

        function redraw!()
            t_start = time(); tmark[] = t_start
            redrawing[] = true
            try
                Z      = ELEMENTS.Z[findfirst(==(elem_dd.value[]), OPTS)]
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

                key = (Z, Tuple(ions), medium, source, T, saha, saha ? ne_exp : 0.0)
                if cache[] === nothing || cache[].key != key
                    spec = prepare_lines(DF, Z; ions = ions, medium = medium, source = source, T = T,
                                         saha = saha, ne = 10.0^ne_exp)
                    cache[] = (key = key, spec = spec)
                end
                lines   = cache[].spec.lines
                use_lte = cache[].spec.use_lte
                lap("odczyt ustawień + linie (cache)")

                # skala osi X (przełączana, gdy na osiach nie ma wykresów)
                if logmode != log_x[]
                    drop_segments!()    # zmieniamy skalę, gdy na osiach nie ma wykresów
                    if logmode
                        # Makie waliduje obecne limity przy zmianie skali, a (0, 3000) nie nadaje się
                        # do log10, więc przed przełączeniem ustawiamy dodatnie
                        lo_s = max(view_x[][1], X_MIN_LOG)
                        hi_s = max(view_x[][2], lo_s * 1.001)
                        xlims!(ax, lo_s, hi_s)
                        xlims!(ax_strip, lo_s, hi_s)
                    end
                    ax.xscale = logmode ? log10 : identity
                    ax_strip.xscale = logmode ? log10 : identity
                    log_x[] = logmode      # dopiero po udanej zmianie
                end
                x0, x1 = clamp_view(view_x[]..., logmode)
                view_x[] = (x0, x1)

                res = max(Float64(res_sl.value[]), RES_FLOOR)
                if rescb.value[]
                    # rozdzielczość = px_sl pikseli ekranu przy bieżącym powiększeniu
                    wpx = max(Float64(ax.scene.viewport[].widths[1]), 200.0)
                    px  = Float64(px_sl.value[])
                    res = max(logmode ? px * log(x1 / x0) / wpx * sqrt(x0 * x1) :
                                        px * (x1 - x0) / wpx, RES_FLOOR)
                end
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
                    (cr, cg, cb) = bw ? (1.0, 1.0, 1.0) : wl2rgb(cl.λ[k])
                    a = 0.15 + 0.85 * h
                    i2 = 2 * j
                    pts1[i2 - 1] = Point2f(cl.λ[k], 0);  pts1[i2] = Point2f(cl.λ[k], 1)
                    col1[i2 - 1] = RGBAf(cr, cg, cb, a);  col1[i2] = RGBAf(cr, cg, cb, a)
                    pts2[i2 - 1] = Point2f(cl.λ[k], 0);  pts2[i2] = Point2f(cl.λ[k], h)
                    col2[i2 - 1] = RGBAf(cr, cg, cb, 1);  col2[i2] = RGBAf(cr, cg, cb, 1)
                end
                lap("punkty i kolory ($nseg kresek)")
                sg.pts1[] = pts1;  sg.col1[] = col1
                sg.pts2[] = pts2;  sg.col2[] = col2
                lap("aktualizacja danych wykresów")
                xlabel_txt = "λ [nm, $(medium == "air" ? "powietrze" : "próżnia")]"
                ax.xlabel[] == xlabel_txt || (ax.xlabel = xlabel_txt)
                xlims!(ax, x0, x1)
                ax.xticks = nice_ticks(x0, x1, logmode)
                ylims!(ax, 0, 1.05)
                ylims!(ax_strip, 0, 1)
                lap("limity i ticki")

                shown[] = (λ = cl.λ[idx], I = cl.I[idx], Iref = Iref, rows = cl.rows[idx], n = cl.n[idx], spread = cl.spread[idx], lines = lines)
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
        onany(elem_dd.value, med_dd.value, xs_dd.value, scale_dd.value, src_dd.value, ref_dd.value,
              ion1_cb.value, ion2_cb.value, ion3_cb.value, ion4_cb.value, ion5_cb.value,
              saha_cb.value, ne_sl.value, bw_cb.value, T_sl.value, res_sl.value, thr_sl.value, rescb.value, px_sl.value) do args...
            request_redraw!()
        end

        # powiększanie/pomniejszanie (factor) i przesuwanie (shift, ułamek nowej szerokości okna);
        # w skali logarytmicznej liczymy w log10, więc zachowujemy stosunek hi/lo
        function move_view!(factor, shift)
            logm = log_x[]
            t  = logm ? log10 : identity
            ti = logm ? (x -> 10.0^x) : identity
            lo, hi = view_x[]
            a0 = t(max(lo, logm ? X_MIN_LOG : X_MIN))
            b0 = t(hi)
            w  = (b0 - a0) * factor
            a  = (a0 + b0) / 2 - w / 2 + shift * w
            view_x[] = clamp_view(ti(a), ti(a + w), logm)
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
                Bonito.DOM.p("Obsługa: kółko myszy = zoom, przeciągnięcie lewym = zaznaczenie obszaru, " *
                             "prawym = przesuwanie; przyciski ◀ ▶ ＋ － działają bez myszy.");
                style = "max-width:760px;margin:8px 0;padding:10px 14px;background:#222228;" *
                        "border:1px solid #444;border-radius:6px;font-size:0.92em;line-height:1.4;")
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
            lo, hi = limits_to_nm(r)
            clo, chi = clamp_view(lo, hi, log_x[])
            if abs(clo - lo) > 1e-9 * max(1.0, abs(lo)) || abs(chi - hi) > 1e-9 * max(1.0, abs(hi))
                redrawing[] = true
                try
                    xlims!(ax, clo, chi)
                finally
                    redrawing[] = false
                end
            end
            view_x[] = (clo, chi)
            ax.xticks = nice_ticks(clo, chi, log_x[])
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
                return
            end
            mraw = mouseposition(ax.scene)[1]
            mx = (log_x[] && mraw < 4.0) ? 10.0^mraw : mraw    # j.w.: log10 -> nm
            if log_x[] && mx <= 0
                return
            end
            t = log_x[] ? log10 : identity
            x0, x1 = view_x[]
            tol = 6 * (t(x1) - t(x0)) / ax.scene.viewport[].widths[1]      # ok. 6 pikseli
            k = searchsortedfirst(λs, mx)
            cand = [i for i in (k - 1, k) if 1 <= i <= length(λs)]
            i = cand[argmin([abs(t(λs[j]) - t(mx)) for j in cand])]
            if abs(t(λs[i]) - t(mx)) <= tol
                set_hover!(describe_line(s.lines, s.rows[i]; n = s.n[i], spread = s.spread[i]))
                set_summary!("$(lines_pl(s.n[i])) o względnej sumie natężeń $(fmt_rel(100 * s.I[i] / s.Iref))",
                             s.n[i] > 1 ? "Możesz zaznaczyć i powiększyć ten wycinek w celu zwiększenia rozdzielczości." : "")
            else
                set_hover!("Najedź kursorem na linię.")
                set_summary!("", "")
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
        drag_x0  = Ref(1.0)
        function mouse_nm()
            mraw = mouseposition(ax.scene)[1]
            return (log_x[] && mraw < 4.0) ? 10.0^mraw : mraw      # j.w.: log10 -> nm
        end
        in_axes() = Makie.is_mouseinside(ax.scene) || Makie.is_mouseinside(ax_strip.scene)
        # wysoki priorytet: Makie przechwytuje zdarzenia myszy w trakcie przeciągania (własny zoom),
        # więc listenery o priorytecie 0 nie dostawałyby ruchów ani puszczenia przycisku
        on(events(fig).mousebutton; priority = 10_000) do ev
            ev.button == Makie.Mouse.left || return
            if ev.action == Makie.Mouse.press
                if in_axes()
                    x = mouse_nm()
                    if x > 0
                        dragging[] = true
                        drag_x0[] = x
                        sel[] = (x, x)
                    end
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
            x = mouse_nm()
            x > 0 || return
            lo, hi = view_x[]
            sel[] = (drag_x0[], clamp(x, lo, hi))
            sel_vis[] || (sel_vis[] = true)
        end

        redraw!()    # pierwsze rysowanie

        # ---------------- układ strony ----------------
        row(items...) = Bonito.DOM.div(items...;
            style = "display:flex;flex-wrap:wrap;align-items:center;gap:16px;margin-bottom:8px;")
        lab(text, items...) = Bonito.DOM.span(Bonito.DOM.span(text * " "), items...;
            style = "display:inline-flex;align-items:center;gap:6px;")
        val(obs) = Bonito.DOM.span(map(string, obs); style = "min-width:3em;")

        return Bonito.DOM.div(
            Bonito.DOM.h2("Spektruś.jl"; style = "margin:0 0 10px 0;"),
            row(lab("Pierwiastek:", elem_dd),
                lab("Zakres:", btns...),
                lab("Widok:", b_left, b_right, b_zin, b_zout),
                lab("Długość fali:", med_dd),
                lab("Skala λ:", xs_dd)),
            row(Bonito.DOM.span("Rysuj linie:"),
                lab("niezjonizowane (I)", ion1_cb),
                lab("zjonizowane (II)", ion2_cb),
                lab("podwójnie (III)", ion3_cb),
                lab("potrójnie (IV)", ion4_cb),
                lab("poczwórnie (V)", ion5_cb),
                lab("czarno-białe", bw_cb),
                about_btn),
            Bonito.DOM.div(about_panel),
            row(lab("Źródło natężeń:", src_dd),
                lab("T [K]:", T_sl, val(T_sl.value)),
                lab("Saha, lg n_e [cm⁻³]:", saha_cb, ne_sl, val(ne_sl.value)),
                lab("Skala natężeń:", scale_dd),
                lab("Odniesienie:", ref_dd),
                lab("Próg [dekady]:", thr_sl, val(thr_sl.value)),
                lab("Rozdzielczość [nm]:", res_sl, val(res_sl.value)),
                lab("lub zależna od zoomu [px]:", rescb, px_sl, val(px_sl.value))),
            Bonito.DOM.div(fig),
            Bonito.DOM.div(info; style = "margin-top:6px;color:#9aa;"),
            Bonito.DOM.div(
                Bonito.DOM.div("Opis linii spektralnych"; style = "font-weight:600;margin-bottom:4px;color:#9ab;"),
                Bonito.DOM.div(summary; style = "min-height:1.4em;"),
                Bonito.DOM.div(hover; style = "font-family:monospace;color:#bbb;min-height:1.4em;margin-top:2px;"),
                Bonito.DOM.div(hint; style = "color:#889;font-size:0.9em;min-height:1.2em;margin-top:2px;");
                style = "margin-top:8px;padding:8px 12px;border:1px solid #444;border-radius:6px;");
            style = "background:#161618;color:#ddd;font-family:system-ui,sans-serif;padding:14px;min-height:100vh;")
    end
end

app = build_app()

# get_server() sam wykrywa port (ENV["PORT"] lub 8080), binduje na 0.0.0.0
# i konfiguruje proxy_url, gdy działa za reverse proxy (JuliaHub, Heroku itp.).
server = Bonito.get_server()
route!(server, "/" => app)
println("Spektruś działa: ", Bonito.online_url(server, "/"))

# Proces trzymamy przy życiu – serwer HTTP działa w tle.
try
    while true
        sleep(1)
    end
catch e
    e isa InterruptException || rethrow()
end