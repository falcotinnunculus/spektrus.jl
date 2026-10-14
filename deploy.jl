using Bonito, WGLMakie

try
    # 1. Tworzymy fizyczny katalog o nazwie "public"
    mkpath("public")
    
    # 2. Wczytujemy aplikację
    include("app.jl")
    println("Zasoby wczytane pomyślnie.")
    
    # 3. Eksportujemy bezpośrednio do pliku index.html wewnątrz folderu public
    export_static("public/index.html", app)
    println("Wyeksportowano aplikację do public/index.html")
catch e
    @error "Błąd podczas eksportu aplikacji!" exception=(e, catch_backtrace())
    rethrow(e)
end