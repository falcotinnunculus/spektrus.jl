using Bonito, WGLMakie

try
    # 1. Tworzymy fizyczny katalog "public"
    mkpath("public")
    
    # 2. Wczytujemy aplikację z app.jl (która tworzy zmienną `app`)
    include("app.jl")
    println("Zasoby wczytane pomyślnie.")
    
    # 3. Eksportujemy statycznie do pliku index.html wewnątrz katalogu public
    export_static("public/index.html", app)
    println("Wyeksportowano aplikację do public/index.html")
catch e
    @error "Błąd podczas eksportu aplikacji!" exception=(e, catch_backtrace())
    rethrow(e)
end