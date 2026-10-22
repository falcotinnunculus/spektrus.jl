# build_sysimage.jl - buduje własny obraz systemowy Julii (sysimage) z wkompilowanymi pakietami aplikacji,
# żeby start i pierwsza strona nie płaciły za kompilację (JIT) Bonito, WGLMakie, CSV i DataFrames.
#
# Wymagania: PackageCompiler zainstalowany w środowisku globalnym (nie w projekcie)
#   julia -e 'using Pkg; Pkg.add("PackageCompiler")'
# oraz kompilator C (gcc lub clang; w obrazie Dockera jest). Na Windowsie to bywa kłopotliwe, dlatego
# sysimage budujemy w Dockerze (patrz Dockerfile i deploy/README.md).
#
# Uruchomienie z folderu projektu (musi tam być spectrum.csv):
#   julia --project=. deploy/build_sysimage.jl
# Potem aplikacja:
#   julia --project=. -J spektrus_sysimage.so app.jl      (.dll na Windowsie, .dylib na macOS)

using PackageCompiler, Pkg, Libdl

root   = normpath(joinpath(@__DIR__, ".."))
wanted = ["Bonito", "WGLMakie", "CSV", "DataFrames", "HTTP"]
deps   = keys(Pkg.project().dependencies)
pkgs   = filter(in(deps), wanted)
missing_pkgs = setdiff(wanted, pkgs)
isempty(missing_pkgs) || @warn "Tych pakietów nie ma w Project.toml, pomijam" missing_pkgs

path = joinpath(root, "spektrus_sysimage." * Libdl.dlext)
@info "Buduję sysimage (to potrwa kilka do kilkunastu minut)" pkgs path
t0 = time()
create_sysimage(pkgs;
    sysimage_path = path,
    precompile_execution_file = joinpath(@__DIR__, "precompile_app.jl"),
    # cpu_target = PackageCompiler.default_app_cpu_target())      # przenośny między procesorami (obraz Dockera działa na różnych maszynach)
    cpu_target="native",
    filter_stdlibs=true)
@info "Gotowe" sekundy = round(Int, time() - t0) rozmiar_MB = round(filesize(path) / 2^20; digits = 1)
println("Uruchomienie: julia --project=. -J $(basename(path)) app.jl")
