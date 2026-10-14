using Bonito, WGLMakie

# Wczytujemy plik app.jl. 
# Ponieważ uruchamiamy deploy.jl, zmienna PROGRAM_FILE będzie inna niż app.jl,
# więc blok z lokalnym serwerem (while true) w app.jl zostanie pominięty!
include("app.jl")

# Wywołujemy eksport statyczny do folderu public
export_static("public", app)
println("Wyeksportowano aplikację do folderu public/")