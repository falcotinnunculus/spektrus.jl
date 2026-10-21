# precompile_app.jl - uruchamiany przez PackageCompiler podczas budowania sysimage (patrz build_sysimage.jl).
# Kompilacje wywołane tym kodem trafiają do obrazu, więc wykonujemy tu to, co aplikacja robi przy starcie i przy
# otwarciu strony: wczytanie danych, rdzeń, budowę strony (handler Bonito) i jej wyrenderowanie.
# Interakcje przez WebSocket (zoom, klikanie) wymagają przeglądarki, więc ich tu nie ma.
#
# Uwaga: do obrazu trafiają skompilowane metody z pakietów (Bonito, Makie, DataFrames itd.). Własny kod aplikacji
# (app.jl, Spektrus.jl) nie jest pakietem, więc jego kompilacja dalej dzieje się przy starcie (robi to warm_core!).

using Bonito, WGLMakie, HTTP

include(joinpath(@__DIR__, "..", "app.jl"))      # wczytuje dane i buduje aplikację; serwer się nie uruchamia (to nie jest plik główny)

warm_core!()

port = 8200 + rand(1:700)
server = Bonito.Server(app, "127.0.0.1", port)
try
    for target in ("/", "/?z=26&ions=1,2&u=4&xs=2&lang=pl", "/?z=11&lang=en")
        HTTP.get("http://127.0.0.1:$port$target"; readtimeout = 300, retry = false)
    end
finally
    close(server)
end
