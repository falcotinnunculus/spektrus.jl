# Spektruś.jl

Widma emisyjne pierwiastków 1–99 w przeglądarce (Bonito + WGLMakie), odtworzenie programu Spektruś 1.0.
Dane linii: NIST Atomic Spectra Database (stopnie jonizacji I–V, 100–2500 nm).

## Pliki

| plik | do czego |
|---|---|
| `app.jl` | aplikacja w przeglądarce |
| `Spektrus.jl` | rdzeń: wczytanie danych, model LTE + Saha, klastrowanie linii, jednostki, nazwy pierwiastków |
| `fetch_nist_lines.jl` | pobiera linie z NIST do katalogu `raw/` (cache) i składa `lines.csv` |
| `prepare_lines.jl` | `lines.csv` → `spectrum.csv` (plik, który czyta aplikacja) |
| `ionization_energies.csv` | energie jonizacji stopni I–V (cm⁻¹, eV), potrzebne do równowagi Saha |
| `test_core.jl` | testy rdzenia |

`spectrum.csv` i `lines.csv` są duże (dziesiątki MB), więc nie ma ich w repozytorium.

## Kolejność uruchamiania

Z folderu projektu (PowerShell):

```
$env:JULIA_NUM_PRECOMPILE_TASKS = 1          # tylko przy pierwszej instalacji pakietów
julia --project=. .\fetch_nist_lines.jl       # raz: pobiera dane z NIST (długo; wynik w raw/ i lines.csv)
julia --project=. .\prepare_lines.jl          # raz po pobraniu: tworzy spectrum.csv
julia --project=. .\test_core.jl              # opcjonalnie: testy rdzenia
julia --project=. .\app.jl                    # aplikacja: http://127.0.0.1:8080  (Ctrl+C kończy)
```

Pobieranie i `prepare_lines.jl` wystarczy powtórzyć tylko wtedy, gdy zmienia się zakres danych (np. stopnie jonizacji).

## Ustawienia w `app.jl`

- granice osi i najmniejsze okno zoomu: `X_MIN`, `X_MAX`, `X_MIN_WIDTH` (stałe na początku pliku),
- najmniejsza rozdzielczość zlewania linii: `RES_FLOOR`,
- wypisywanie czasów przerysowania w konsoli: `TIMING = true` (domyślnie wyłączone).
- język interfejsu na starcie: `LANG_DEFAULT` (`"pl"` albo `"en"`); w aplikacji przełącza go przycisk PL | EN,
- imię i nazwisko autora w oknie „O programie”: `AUTHOR`.

## Jednostki osi

nm, Å, µm, cm⁻¹, eV. Dla cm⁻¹ i eV używane są długości fal w próżni, a oś rośnie z energią (widmo jest odbite względem osi w nm).

## Model natężeń

Emisja optycznie cienka w LTE, `I ∝ (N_z/N)·g_k·A_ki·ν·exp(−E_k/kT)/U_z(T)`.
Udziały stopni `N_z/N` z równowagi Saha dla zadanej gęstości elektronowej (suwak lg n_e), uwzględnione stopnie I–V.
Funkcje podziału `U_z` liczone są z poziomów występujących w pobranych liniach, więc to przybliżenie.

## Cytowanie danych

A. Kramida, Yu. Ralchenko, J. Reader i NIST ASD Team, *NIST Atomic Spectra Database*, https://physics.nist.gov/asd