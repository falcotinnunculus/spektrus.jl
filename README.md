# Spektruś.jl

Widma emisyjne pierwiastków 1–99 w przeglądarce (Bonito + WGLMakie). Odtworzenie programu Spektruś 1.0
(Mikołaj Pytel, 2002) z danymi z NIST Atomic Spectra Database.

Repozytorium: https://github.com/falcotinnunculus/spektrus.jl

## Funkcje

- Linie emisyjne pierwiastków 1–99, stopnie jonizacji I–V, 100–2500 nm (UV, światło widzialne, podczerwień).
- Natężenia z modelu LTE i równowagi Sahy (temperatura T, gęstość elektronowa n_e) albo umowne natężenia NIST.
- Zoom (zaznaczenie, kółko myszy, przyciski), przesuwanie, jednostki nm, Å, µm, cm⁻¹ i eV, skala osi liniowa i logarytmiczna.
- Opis linii po najechaniu (przejście, Aki, natężenie), odczyt położenia kursora, tabela linii z widoku z eksportem CSV i TSV.
- Zapis wykresu jako PNG (ciemne tło albo jasne, do druku) i link do widoku (ustawienia są w adresie strony).
- Gotowe przykłady (lampy wodorowa, helowa, neonowa, sodowa, rtęciowa, wapń z liniami H i K), pomoc ze wstępem do spektroskopii, skróty klawiszowe.
- Interfejs po polsku i angielsku (język z przeglądarki).

## Wymagania i instalacja

Julia 1.11 (testowane na 1.11.4). Pakiety: Bonito, WGLMakie, CSV, DataFrames oraz HTTP (tylko skrypty pobierające dane).
Z folderu projektu (PowerShell):

```
$env:JULIA_NUM_PRECOMPILE_TASKS = 1          # tylko przy pierwszej instalacji pakietów (mniej RAM-u)
julia --project=.
]instantiate                                  # gdy jest Project.toml; inaczej: ]add Bonito WGLMakie CSV DataFrames HTTP
```

## Dane i kolejność uruchamiania

| plik | do czego |
|---|---|
| `app.jl` | aplikacja w przeglądarce |
| `Spektrus.jl` | rdzeń: dane, model LTE + Saha, łączenie linii, jednostki, nazwy pierwiastków, tabela linii |
| `fetch_nist_lines.jl` | pobiera linie z NIST do `raw/` (cache) i składa `lines.csv` |
| `prepare_lines.jl` | `lines.csv` → `spectrum.csv` (plik, który czyta aplikacja) |
| `fetch_nist_levels.jl` | (opcjonalnie) pobiera pełne listy poziomów z NIST do `raw_levels/` i składa `levels.csv` (dokładniejsze funkcje podziału) |
| `compare_levels.jl` | (narzędzie) porównuje funkcje podziału i udziały stopni z `levels.csv` i bez niego |
| `ionization_energies.csv` | energie jonizacji stopni I–V (cm⁻¹, eV), potrzebne do równowagi Saha |
| `test_core.jl` | testy rdzenia |
| `NOTES.md` | notatki dla programisty: architektura, decyzje i pułapki |

Pliki generowane (`spectrum.csv`, `lines.csv`, `levels.csv`, `raw/`, `raw_levels/`) są duże albo odtwarzalne, więc nie ma ich w repozytorium.

```
julia --project=. .\fetch_nist_lines.jl       # raz: pobiera dane z NIST (długo; wynik w raw/ i lines.csv)
julia --project=. .\prepare_lines.jl          # raz po pobraniu: tworzy spectrum.csv
julia --project=. .\fetch_nist_levels.jl      # opcjonalnie: levels.csv (ok. 10 min; zacznij od `... 26 26` dla samego żelaza)
julia --project=. .\test_core.jl              # opcjonalnie: testy rdzenia
julia --project=. .\app.jl                    # aplikacja: http://127.0.0.1:8080  (Ctrl+C kończy)
```

Pobieranie i `prepare_lines.jl` wystarczy powtórzyć tylko wtedy, gdy zmienia się zakres danych (np. stopnie jonizacji).
Przy starcie aplikacja przez kilkanaście sekund przygotowuje rdzeń (kompilacja, wczytanie `levels.csv`), dzięki czemu pierwsza strona
ładuje się szybciej.

## Obsługa

- **Wykres:** kółko myszy = zoom, przeciągnięcie lewym przyciskiem = zaznaczenie obszaru do powiększenia, prawym = przesuwanie.
- **Skróty klawiszowe:** ← → przesuwanie, + − zoom, Home światło widzialne, End całe widmo, F1 pomoc, Esc zamyka okna.
  Nie działają, gdy aktywne jest pole tekstowe, lista lub suwak.
- **Pomoc:** przycisk „Pomoc” (F1) zawiera instrukcję, krótkie wprowadzenie do spektroskopii i opis modelu natężeń.
- **Link do widoku:** adres strony zawiera ustawienia, więc odświeżenie strony je zachowuje, a „Link” kopiuje adres.

Parametry adresu (wszystkie opcjonalne, np. `/?z=26&ions=1,2&x=380,780&T=8000`):

| parametr | znaczenie |
|---|---|
| `z` | liczba atomowa (1–99) |
| `x` | zakres osi `od,do` w nm |
| `u` | jednostka: 1 nm, 2 Å, 3 µm, 4 cm⁻¹, 5 eV |
| `xs` | skala osi: 1 liniowa, 2 logarytmiczna |
| `sc` | skala natężeń: 1 logarytmiczna, 2 liniowa |
| `src` | źródło natężeń: 1 model LTE, 2 NIST |
| `ref` | odniesienie jasności: 1 linia widzialna, 2 linia w oknie, 3 linia pierwiastka |
| `med` | długość fali w: 1 powietrzu, 2 próżni |
| `ions` | stopnie jonizacji, np. `1,2` |
| `T`, `ne` | temperatura [K] i lg n_e [cm⁻³] |
| `thr` | zakres jasności [dekady] |
| `res`, `px`, `rescb` | rozdzielczość zlewania [nm], rozdzielczość w pikselach, `rescb=1`: zależna od zoomu |
| `saha`, `bw`, `gray`, `mem` | `1`/`0`: równowaga Sahy, czarno-białe, szare linie poza światłem widzialnym, składowe zlanych linii w opisie |
| `lang` | `pl` albo `en` (ma pierwszeństwo przed językiem przeglądarki) |

## Ustawienia w kodzie (`app.jl`)

| stała | znaczenie |
|---|---|
| `X_MIN`, `X_MIN_LOG`, `X_MAX`, `X_MIN_WIDTH` | granice osi i najmniejsze okno zoomu |
| `RES_FLOOR` | najmniejsza rozdzielczość zlewania linii [nm] |
| `TABLE_ROWS` | ile linii pokazuje tabela (eksport zawiera wszystkie) |
| `LANG_DEFAULT` | język, gdy przeglądarka nie podaje polskiego ani angielskiego (domyślnie `"en"`) |
| `AUTHOR` | autor w oknie „O programie” |
| `TIMING`, `LOAD_TIMING` | czasy przerysowania wykresu i budowania strony w konsoli (domyślnie wyłączone) |
| `PRESETS` (w `build_app`) | przykłady: pierwiastek, temperatura i stopnie jonizacji; temperatury są poglądowe |

Zmienna środowiskowa `SPEKTRUS_WARM=0` wyłącza przygotowanie rdzenia przy starcie (PowerShell: `$env:SPEKTRUS_WARM = "0"`).

## Jednostki osi

nm, Å, µm, cm⁻¹, eV. Dla cm⁻¹ i eV używane są długości fal w próżni, a oś rośnie z energią (widmo jest odbite względem osi w nm).

## Model natężeń

Emisja optycznie cienka w LTE, `I ∝ (N_z/N)·g_k·A_ki·ν·exp(−E_k/kT)/U_z(T)`.
Udziały stopni `N_z/N` z równowagi Saha dla zadanej gęstości elektronowej (suwak lg n_e), uwzględnione stopnie I–V.
Funkcje podziału `U_z` liczone są z pełnych list poziomów NIST (`levels.csv`) uzupełnionych poziomami z pobranych linii
(część jonów ma w NIST tylko poziom podstawowy). Gdy `levels.csv` nie ma, używane są tylko poziomy występujące w pobranych
liniach, co jest przybliżeniem (dla złożonych atomów może zaniżać `U`).
Linie bez współczynnika Aki nie mają natężenia w modelu LTE: gdy żaden z wybranych stopni nie ma Aki, program używa
umownych natężeń NIST, a gdy tylko część, pomija te stopnie i ostrzega o tym pod wykresem.
To model poglądowy (nauka, porównania), nie narzędzie pomiarowe.

## Testy

```
julia --project=. .\test_core.jl
```

Testy sprawdzają rdzeń (`Spektrus.jl`) na prawdziwych danych. Interfejs (`app.jl`) nie ma testów automatycznych i jest sprawdzany ręcznie.

## Źródła i cytowanie

- Dane linii i poziomów: A. Kramida, Yu. Ralchenko, J. Reader i NIST ASD Team, *NIST Atomic Spectra Database*, https://physics.nist.gov/asd
- M. N. Saha, *Ionization in the solar chromosphere*, Philosophical Magazine 40, 472 (1920), doi:10.1080/14786441008636148
- H. R. Griem, *Principles of Plasma Spectroscopy*, Cambridge University Press, 1997 (rozdz. 7: Thermodynamic equilibrium relations)
- Kolory linii: przybliżenie długości fali na RGB według algorytmu Dana Brutona (jak w oryginalnym Spektrusiu).

## Autorzy

Mateusz Winiarski i Claude (Anthropic).