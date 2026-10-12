using DataFrames
using Gtk4

# Utwórz główne okno
window = GtkWindow("Spektrus.jl", 600, 400)

# Utwórz kontener pionowy
vbox = GtkBox(:v, 5)
set_gtk_property!(vbox, :margin, 10)

# Utwórz przycisk i etykietę
button = GtkButton("Kliknij mnie!")
label = GtkLabel("Naciśnij przycisk")

# Dodaj elementy do kontenera
push!(vbox, button)
push!(vbox, label)
push!(window, vbox)

# Funkcja obsługi kliknięcia przycisku
function button_clicked_callback(widget)
    set_gtk_property!(label, :label, "Witaj Świecie!")
end

# Połączenie sygnału z funkcją callback (nowa składnia)
signal_connect(button_clicked_callback, button, "clicked")

# Wyświetl okno (zaktualizowane z showall)
show(window)

# Uruchom główną pętlę GTK
if !isinteractive()
    Gtk.main()
end