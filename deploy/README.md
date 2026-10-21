# Wdrożenie

Aplikacja to serwer Julii z połączeniem WebSocket na każdą otwartą kartę, więc potrzebuje hosta, na którym działa proces
(nie nadaje się na statyczny hosting). Poniżej droga przez obraz Dockera ze sysimage: ten sam obraz działa lokalnie,
na Cloud Run i na zwykłym VPS.

## 1. Dane jako plik w GitHub Releases

`spectrum.csv` i `levels.csv` nie są w repozytorium, a pobranie z NIST trwa długo, więc obraz pobiera je z gotowego pliku.
W folderze z danymi (PowerShell ma wbudowane `tar`):

```
tar -czf spektrus-data.tar.gz spectrum.csv levels.csv
```

Na GitHubie: Releases → Draft a new release → tag `data-v1` → dołącz `spektrus-data.tar.gz` → Publish. Adres pliku to
`https://github.com/falcotinnunculus/spektrus.jl/releases/download/data-v1/spektrus-data.tar.gz`
(taki jest domyślny `DATA_URL` w `Dockerfile`; przy innym tagu: `docker build --build-arg DATA_URL=...`).
Dane pochodzą z NIST ASD, więc zostaw cytowanie z README.

## 2. Test lokalny w Dockerze

Potrzebny `Project.toml` i `Manifest.toml` w repozytorium (`Manifest.toml` nie jest już w `.gitignore`: wdrożenie ma budować
dokładnie te wersje pakietów, które testowałeś).

```
docker build -t spektrus .                 # kilka do kilkunastu minut (to szacunek): głównie budowanie sysimage
docker run --rm -p 8080:8080 spektrus      # http://localhost:8080
```

Sprawdź w konsoli „Przygotowanie rdzenia: X s” i czas pierwszej strony, bo dopiero to pokaże, ile daje sysimage. Pierwsze
budowanie może się nie powieść (sysimage jest wrażliwy na wersje pakietów): w takim razie wklej mi koniec logu.

## 3. Google Cloud Run (darmowy przy małym ruchu)

Wymaga konta Google Cloud z włączonymi płatnościami (karta), choć przy małym ruchu rachunek ma wynosić 0. Darmowy limit
miesięczny (rozliczanie za żądania, wg cennika Google): 180 000 vCPU-sekund, 360 000 GiB-sekund i 2 mln żądań. Przy 1 vCPU
i 2 GiB pamięci to ok. 50 godzin połączonych sesji miesięcznie. Połączenie WebSocket liczy się jako trwające żądanie,
więc otwarta karta zużywa limit. Po jego przekroczeniu zaczyna się naliczanie opłat.

```
gcloud auth login
gcloud config set project TWOJ_PROJEKT
gcloud services enable run.googleapis.com cloudbuild.googleapis.com artifactregistry.googleapis.com

gcloud run deploy spektrus --source . --region europe-west1 \
  --memory 2Gi --cpu 1 --cpu-boost --min-instances 0 --max-instances 1 \
  --concurrency 8 --timeout 3600 --allow-unauthenticated
```

- `--min-instances 0`: brak opłat bez ruchu, ale po przerwie pierwsza osoba czeka na zimny start (kontener, start Julii,
  wczytanie danych, przygotowanie rdzenia). Ile to trwa w Twoim przypadku, pokaże dopiero pomiar.
- `--max-instances 1` i `--concurrency 8`: ograniczają koszt. Gdy ktoś zostawi otwarte karty, instancja nie rozrasta się.
- `--timeout 3600`: maksymalny czas połączenia (godzina); po nim przeglądarka musi się połączyć ponownie.
- Region: wybierz region z niższym cennikiem (Tier 1, np. `europe-west1`); darmowy limit jest równoważnikiem kwotowym,
  więc w droższym regionie starczy na mniej godzin. Sprawdź listę w cenniku Cloud Run.
- Jeśli budowanie z `--source` przekroczy domyślny limit czasu Cloud Build (domyślnie ok. 10 minut, sprawdź), zbuduj
  obraz osobno: `gcloud builds submit --timeout=3600s --tag REGION-docker.pkg.dev/PROJEKT/REPO/spektrus`
  i wdróż go przez `gcloud run deploy spektrus --image ...` z tymi samymi parametrami.

Po pierwszym wdrożeniu Cloud Run poda adres usługi (`https://spektrus-XXXX.europe-west1.run.app`). Bonito musi go znać,
żeby przeglądarka łączyła się WebSocketem we właściwe miejsce:

```
gcloud run services update spektrus --region europe-west1 \
  --update-env-vars SPEKTRUS_URL=https://spektrus-XXXX.europe-west1.run.app
```

**Ustaw alert budżetowy** (Console → Billing → Budgets & alerts, np. 1 USD), żeby zauważyć, gdy limit darmowy się skończy.

Jeśli strona się ładuje, ale wykres się nie pojawia albo połączenie ciągle się zrywa, to najczęściej zły `SPEKTRUS_URL`
(adres musi zaczynać się od `https://` i nie mieć ukośnika na końcu).

## 4. Gdy zimny start przeszkadza: VPS

Najtańszy VPS z 2 GB RAM kosztuje kilka euro miesięcznie (nie jest darmowy) i działa cały czas, bez zimnych startów.
Na serwerze z Dockerem i Caddy (sam załatwia HTTPS i obsługuje WebSocket):

```
docker run -d --restart unless-stopped --name spektrus -p 127.0.0.1:8080:8080 \
  -e SPEKTRUS_URL=https://spektrus.twojadomena.pl spektrus
```

`Caddyfile`:

```
spektrus.twojadomena.pl {
    reverse_proxy 127.0.0.1:8080
}
```

## 5. Czego nie sprawdziłem

Nie uruchamiałem ani budowania obrazu, ani sysimage, ani aplikacji za serwerem pośredniczącym (nie mam tu Julii ani Dockera).
Nie znam zużycia pamięci procesu ani czasu zimnego startu. Sysimage obejmuje skompilowane metody pakietów (Bonito, Makie,
DataFrames itd.), a własny kod aplikacji nadal się kompiluje przy starcie (robi to `warm_core!`), więc spodziewam się dużej,
ale nie całkowitej poprawy. Wątek na forum Julii z kwietnia 2026 opisuje kłopoty z WGLMakie za nginxem, więc test za
serwerem pośredniczącym (Cloud Run, Caddy) jest konieczny.
