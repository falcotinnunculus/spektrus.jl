# bin/main.jl

# 1. Przechodzimy do katalogu głównego projektu.
# Dzięki temu ścieżki do plików typu catalog.dat, lines.csv, spectrum.csv
# będą działać tak samo, jak przy lokalnym uruchamianiu.
cd(joinpath(@__DIR__, ".."))

# 2. Uruchamiamy Twój główny skrypt aplikacji.
include("app.jl")