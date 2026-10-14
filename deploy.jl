using Bonito, WGLMakie

try
    include("app.jl")
    println("Zasoby wczytane pomyślnie.")
    
    export_static("public", app)
    println("Wyeksportowano aplikację do folderu public/")
catch e
    @error "Błąd podczas eksportu aplikacji!" exception=(e, catch_backtrace())
    rethrow(e)
end