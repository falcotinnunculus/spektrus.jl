# Obraz do wdrożenia (Cloud Run, VPS, każdy host z Dockerem). Szczegóły: deploy/README.md
#   docker build -t spektrus .
#   docker run --rm -p 8080:8080 spektrus        (potem http://localhost:8080)
FROM julia:1.11

# gcc jest potrzebny PackageCompiler do zbudowania sysimage
RUN apt-get update \
 && apt-get install -y --no-install-recommends gcc libc6-dev curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# PackageCompiler w środowisku globalnym, żeby nie trafił do Project.toml aplikacji
RUN julia -e 'using Pkg; Pkg.add("PackageCompiler")'

# pakiety aplikacji (osobna warstwa: przebudowuje się tylko po zmianie Project.toml i Manifest.toml)
COPY Project.toml ./
RUN julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()'

COPY app.jl Spektrus.jl ionization_energies.csv ./
COPY deploy ./deploy

# dane (spectrum.csv i levels.csv) jako plik w GitHub Releases, bo nie ma ich w repozytorium
ARG DATA_URL=https://github.com/falcotinnunculus/spektrus.jl/releases/download/data-v1/spektrus-data.tar.gz
RUN curl -fsSL "$DATA_URL" | tar -xz

# sysimage: najdłuższy krok budowania (kilka do kilkunastu minut)
RUN JULIA_IMAGE_THREADS=1 julia --project=. deploy/build_sysimage.jl

ENV HOST=0.0.0.0 PORT=8080
EXPOSE 8080
CMD ["julia", "--project=.", "-J", "spektrus_sysimage.so", "app.jl"]
