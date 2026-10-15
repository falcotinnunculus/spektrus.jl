using Bonito

# Zdefiniuj swoją aplikację
example_app = App(DOM.div("Witaj z Bonito na JuliaHub!"), title="Moja aplikacja")

# get_server() automatycznie wykryje proxy JuliaHub i uruchomi serwer
server = Bonito.get_server()

# Dodaj trasę główną wskazującą na Twoją aplikację
route!(server, "/" => example_app)

# Serwer działa i jest dostępny pod adresem zwracanym przez online_url
println("Aplikacja dostępna pod: ", online_url(server, "/"))