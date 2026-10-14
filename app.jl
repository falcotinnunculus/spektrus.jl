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
include(joinpath(@__DIR__, "spektrus.jl"))
using .Spektrus

WGLMakie.activate!()

# ---- granice osi λ (zmień tutaj) ----
const X_MIN       = 0.0       # dolna granica [nm], skala liniowa
const X_MIN_LOG   = 10.0      # dolna granica [nm], skala logarytmiczna (log(0) nie istnieje)
const X_MAX       = 3000.0    # górna granica [nm]; dane NIST sięgają 2500 nm
const X_MIN_WIDTH = 0.05      # najwęższe okno [nm] w skali liniowej (w logarytmicznej: 0,002 dekady)

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
    w = clamp(t(hi) - a, logmode ? 0.002 : X_MIN_WIDTH, amax - amin)
    if a < amin
        a = amin
    end
    if a + w > amax
        a = amax - w
    end
    return (ti(a), ti(a + w))
end

function build_app()
    return Bonito.App() do session
        # ---------------- kontrolki ----------------
        elem_dd  = Bonito.Dropdown(OPTS; index = something(findfirst(==(11), ELEMENTS.Z), 1))
        med_dd   = Bonito.Dropdown(["powietrze", "próżnia"])
        xs_dd    = Bonito.Dropdown(["liniowa", "logarytmiczna"])
        scale_dd = Bonito.Dropdown(["logarytmiczna", "liniowa"])
        src_dd   = Bonito.Dropdown(["model LTE", "NIST"])
        ref_dd   = Bonito.Dropdown(["linia widzialna", "linia w oknie", "linia pierwiastka"])
        ion_cb   = Bonito.Checkbox(false)
        T_sl     = Bonito.Slider(1000:100:10000; value = 6000)
        res_sl   = Bonito.Slider(0.0:0.1:5.0; value = 0.0)
        thr_sl   = Bonito.Slider(0.5:0.1:6.0; value = 4.0)
        btns     = [Bonito.Button(name) for (name, _) in RANGES]

        info  = Observable("")
        hover = Observable("Najedź kursorem na linię.")

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
        ax = Axis(fig[2, 1]; xlabel = "λ [nm]", ylabel = "natężenie względne", axstyle...)
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

        function redraw!()
            redrawing[] = true
            try
                Z      = ELEMENTS.Z[findfirst(==(elem_dd.value[]), OPTS)]
                medium = med_dd.value[] == "powietrze" ? "air" : "vac"
                source = src_dd.value[] == "model LTE" ? "lte" : "nist"
                scale  = scale_dd.value[] == "logarytmiczna" ? "log" : "lin"
                refmode = ref_dd.value[] == "linia w oknie" ? "window" :
                          ref_dd.value[] == "linia pierwiastka" ? "all" : "visible"
                T       = Float64(T_sl.value[])
                res     = max(Float64(res_sl.value[]), 0.01)
                decades = Float64(thr_sl.value[])
                ion2    = ion_cb.value[]
                logmode = xs_dd.value[] == "logarytmiczna"

                key = (Z, ion2, medium, source, T)
                if cache[] === nothing || cache[].key != key
                    spec = prepare_lines(DF, Z; ion2 = ion2, medium = medium, source = source, T = T)
                    cache[] = (key = key, spec = spec)
                end
                lines   = cache[].spec.lines
                use_lte = cache[].spec.use_lte
                cl      = cluster_window(lines, -Inf, Inf, res)

                empty!(ax)
                empty!(ax_strip)

                # skala osi X (przełączana, gdy na osiach nie ma wykresów)
                if logmode != log_x[]
                    log_x[] = logmode
                    ax.xscale = logmode ? log10 : identity
                    ax_strip.xscale = logmode ? log10 : identity
                end
                x0, x1 = clamp_view(view_x[]..., logmode)
                view_x[] = (x0, x1)

                Iref = reference_intensity(lines, x0, x1; ref = refmode)
                idx  = isempty(cl.I) ? Int[] : findall(visible_mask(cl.I, Iref, decades))

                if !isempty(idx)
                    pts1 = Point2f[]; col1 = RGBAf[]
                    pts2 = Point2f[]; col2 = RGBAf[]
                    for k in idx
                        h = line_height(cl.I[k], Iref, decades; scale = scale)
                        (cr, cg, cb) = wl2rgb(cl.λ[k])
                        a = 0.15 + 0.85 * h
                        push!(pts1, Point2f(cl.λ[k], 0), Point2f(cl.λ[k], 1))
                        push!(col1, RGBAf(cr, cg, cb, a), RGBAf(cr, cg, cb, a))
                        push!(pts2, Point2f(cl.λ[k], 0), Point2f(cl.λ[k], h))
                        push!(col2, RGBAf(cr, cg, cb, 1), RGBAf(cr, cg, cb, 1))
                    end
                    linesegments!(ax_strip, pts1; color = col1, linewidth = 1.5)
                    linesegments!(ax, pts2; color = col2, linewidth = 1.5)
                end
                ax.xlabel = "λ [nm, $(medium == "air" ? "powietrze" : "próżnia")]"
                xlims!(ax, x0, x1)
                ylims!(ax, 0, 1.05)
                ylims!(ax_strip, 0, 1)

                shown[] = (λ = cl.λ[idx], rows = cl.rows[idx], n = cl.n[idx], lines = lines)
                src_txt = use_lte ? "model LTE, T = $(Int(T)) K" :
                          (source == "lte" ? "NIST (brak Aki dla tego pierwiastka, LTE niedostępny)" : "NIST")
                info[] = "Linie widoczne: $(length(idx)) z $(length(cl.λ)) (całe widmo)   ·   " *
                         "źródło natężeń: $src_txt   ·   odniesienie jasności: $(round(Iref; sigdigits = 3))"
            finally
                redrawing[] = false
            end
            return nothing
        end

        # ---------------- reakcje ----------------
        onany(elem_dd.value, med_dd.value, xs_dd.value, scale_dd.value, src_dd.value, ref_dd.value,
              ion_cb.value, T_sl.value, res_sl.value, thr_sl.value) do args...
            redraw!()
        end

        for (b, (_, (lo, hi))) in zip(btns, RANGES)
            on(b.value) do _
                view_x[] = (lo, hi)
                redraw!()
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
            ref_dd.value[] == "linia w oknie" && redraw!()
        end

        # podpowiedź po najechaniu na linię
        on(events(fig).mouseposition) do _
            s = shown[]
            s === nothing && return
            if !(Makie.is_mouseinside(ax.scene) || Makie.is_mouseinside(ax_strip.scene))
                return
            end
            λs = s.λ
            if isempty(λs)
                hover[] = "Brak widocznych linii."
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
                hover[] = describe_line(s.lines, s.rows[i]; n = s.n[i])
            else
                hover[] = "Najedź kursorem na linię."
            end
        end

        redraw!()    # pierwsze rysowanie

        # ---------------- układ strony ----------------
        row(items...) = Bonito.DOM.div(items...;
            style = "display:flex;flex-wrap:wrap;align-items:center;gap:16px;margin-bottom:8px;")
        lab(text, items...) = Bonito.DOM.span(Bonito.DOM.span(text * " "), items...;
            style = "display:inline-flex;align-items:center;gap:6px;")
        val(obs) = Bonito.DOM.span(map(string, obs); style = "min-width:3em;")

        return Bonito.DOM.div(
            Bonito.DOM.h2("Spektruś.jl BETA"; style = "margin:0 0 10px 0;"),
            row(lab("Pierwiastek:", elem_dd),
                lab("Zakres:", btns...),
                lab("Długość fali:", med_dd),
                lab("Skala λ:", xs_dd),
                lab("jony II", ion_cb)),
            row(lab("Źródło natężeń:", src_dd),
                lab("T [K]:", T_sl, val(T_sl.value)),
                lab("Skala natężeń:", scale_dd),
                lab("Odniesienie:", ref_dd),
                lab("Próg [dekady]:", thr_sl, val(thr_sl.value)),
                lab("Rozdzielczość [nm]:", res_sl, val(res_sl.value))),
            Bonito.DOM.div(fig),
            Bonito.DOM.div(info; style = "margin-top:6px;color:#9aa;"),
            Bonito.DOM.div(hover; style = "margin-top:4px;font-family:monospace;min-height:1.5em;");
            style = "background:#161618;color:#ddd;font-family:system-ui,sans-serif;padding:14px;min-height:100vh;")
    end
end

app = build_app()

Bonito.js"""document.title = "spektrus.jl BETA";"""

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
