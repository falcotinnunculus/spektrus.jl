# Notatki dla programisty

Krótki opis budowy programu, decyzji i pułapek, na które natrafiliśmy. Użytkownika interesuje raczej README.

## Przepływ danych

```
fetch_nist_lines.jl  ->  raw/  ->  lines.csv  ->  prepare_lines.jl  ->  spectrum.csv   (czyta app.jl przez Spektrus.load_data)
fetch_nist_levels.jl ->  raw_levels/  ->  levels.csv                                   (opcjonalne; Spektrus.level_db)
ionization_energies.csv                                                                (Spektrus.load_ionization)
```

## Rdzeń (`Spektrus.jl`)

- `prepare_lines(df, Z; ions, medium, source, T, saha, ne)` wybiera linie pierwiastka i liczy natężenia `I`.
  Tryb LTE pomija linie bez Aki, a stopnie, w których pominięto wszystkie linie, zwraca w polu `skipped`.
- `ionization_balance` liczy udziały stopni (Saha) z funkcji podziału `U_z` (`level_table`, `partition_function`).
  `level_table` łączy poziomy z `levels.csv` z poziomami z linii (NIST wygrywa przy tej samej energii).
  Cache `LEVELS` ma klucz `(objectid(df), objectid(baza poziomów), Z, stopień)`.
- `cluster_window` zlewa bliskie linie w jedną kreskę (rozdzielczość `res`), `visible_mask`, `line_height`, `reference_intensity`
  przygotowują wysokości kresek.
- `line_table` buduje tabelę do podglądu i eksportu, `parse_query`, `snap`, `pick_language` obsługują adres i nagłówki.

## Interfejs (`app.jl`)

Cała strona powstaje w jednej funkcji `build_app` (handler `Bonito.App`, osobno dla każdego wczytania strony). Sekcje po kolei:
kontrolki, stan osi i wykres, funkcje rysowania (`redraw!`), reakcje na zmiany, rozgrzewka, tabela, skróty klawiszowe,
zapis PNG, panele (O programie, Pomoc, tabela, menu) i układ.

Rzeczy, które nie są oczywiste:

1. **Napisy PL/EN.** Napisy statyczne są w słowniku `TR` (`tr(klucz)` daje Observable), dynamiczne przez `L(pl, en)`.
   Listy rozwijane Bonito resetują wybór po zmianie opcji, więc każda lista istnieje w dwóch wersjach (`dd2`), ze wspólnym
   indeksem `idx`, a widoczna jest jedna. Kod porównuje numery opcji, nigdy ich tekst.
2. **Szerokości.** `lbl`, `pair` i `dd2` nadają elementom minimalną szerokość (najdłuższy napis z PL i EN, w `em`),
   żeby przełączenie języka nie zmieniało zawijania paska. Przyciski i listy w parze wypełniają swoje miejsce (CSS `.pair`).
3. **Interpolacja w `js"..."`.** Bonito wstawia teksty w apostrofach bez zamiany nowych linii, więc wieloliniowy tekst
   (CSV) przekazujemy jako base64 (`download_text`). Krótkie, jednoliniowe teksty (nazwy plików, adres) można wstawiać wprost.
4. **Zapis PNG.** WebGL czyści bufor po wyświetleniu klatki, więc `export_png!` podpina się pod `requestAnimationFrame` i robi
   zrzut zaraz po narysowaniu klatki. Wariant jasny chwilowo przełącza motyw wykresu (`apply_theme!`) i przyciemnia kolory kresek.
5. **Adres = stan widoku.** `current_qs` składa parametry, `ui_sync!` (z opóźnieniem 0,6 s) wpisuje je do adresu
   (`history.replaceState`) i odświeża tabelę. Startowe wartości czytane są z `request.target` (funkcje `q_num`, `q_int`, `q_bool`, `snap`).
6. **Pole liczby atomowej** to zwykły `<input type=number>` z `Observable{Int}`: `Bonito.NumberInput` istnieje tylko dla `Float64`.
7. **Ostrzeżenia.** Filtr `DropZoomWarn` wycisza nieszkodliwe „Rectangle zoom ignored” (kliknięcie bez przeciągania).
8. **Kompilacja (JIT).** Pierwsze wczytanie strony po starcie Julii jest wolne (kompilacja), kolejne trwają ułamek sekundy.
   `warm_core!` przy starcie skraca to dla rdzenia i Makie, przycisk „Rozgrzej” robi to samo dla interakcji (znika po użyciu).
   Pełne rozwiązanie to sysimage (PackageCompiler), nie jest zrobione.

## Wdrożenie

- `Dockerfile` buduje obraz Julii 1.11, pobiera dane z GitHub Releases i buduje sysimage (`deploy/build_sysimage.jl`, ślad kompilacji
  z `deploy/precompile_app.jl`: wczytanie danych, rdzeń, otwarcie strony przez HTTP). Do sysimage trafiają metody pakietów
  (Bonito, Makie, DataFrames, CSV, HTTP), a własny kod (`app.jl`, `Spektrus.jl`) nie jest pakietem, więc nadal się kompiluje
  przy starcie (`warm_core!`). Przekształcenie `Spektrus` w pakiet pozwoliłoby wkompilować i jego.
- PackageCompiler wymaga kompilatora C (gcc lub clang); na Windowsie bywa kłopotliwy, dlatego sysimage buduje się w Dockerze.
  Sysimage ma `cpu_target` przenośny między procesorami (`default_app_cpu_target`).
- Serwer czyta `HOST`, `PORT` i `SPEKTRUS_URL` (zewnętrzny adres, czyli `proxy_url` w `Bonito.Server`). Za serwerem pośredniczącym
  WebSocket musi być przekazywany; zły adres objawia się pustym wykresem albo zrywanym połączeniem.
- `Manifest.toml` jest w repozytorium, żeby wdrożenie budowało dokładnie przetestowane wersje pakietów.

## Dane

- NIST podaje kilka limitów serii w poziomach jednego jonu; energia jonizacji to najniższy z nich (`fetch_nist_levels.jl`).
- Dla części jonów (np. La III–V, wiele ciężkich pierwiastków) NIST nie podaje Aki, więc w LTE ich linii nie ma.
- Temperatury w `PRESETS` są dobrane poglądowo (nie z literatury): przykłady mają pokazywać charakterystyczne linie.

## Co jest nietestowane automatycznie

Interfejs sprawdzany jest ręcznie w przeglądarce. Lista rzeczy do obejrzenia po większej zmianie: język (PL/EN), przełączanie jednostek,
zoom i przesuwanie (mysz i klawisze), menu i panele, zapis PNG (ciemny i jasny), pobranie CSV i TSV, link do widoku (F5 przywraca stan),
przykłady, zawijanie paska w wąskim oknie.

## Pomysły, których nie zrobiono

- „Co to za linia?”: identyfikacja pierwiastka po długości fali.
- Widmo ciągłe (ciało doskonale czarne) do porównania z liniowym.
- Poprawka na obniżenie energii jonizacji przez plazmę (Debye) przy dużym n_e.
- Oszacowanie Aki dla stopni bez danych (o ile da się z dostępnych danych).
- Sysimage (PackageCompiler) zamiast rozgrzewki.
- Podział `build_app` na mniejsze pliki.
