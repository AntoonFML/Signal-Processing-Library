import os
import platform
import numpy as np
from scipy import signal
import scikit_build_example._core as m

def list_audio_files(directory):
    return [f for f in os.listdir(directory) if f.endswith('.wav')]

def clear_console():
    if platform.system() == "Windows":
        os.system("cls")
    else:
        os.system("clear")

def display_menu():
    print("Wybierz opcje:")
    print("0. Wyswietl wykres sygnalu audio")
    print("1. Wykonaj transformate Fouriera i wyswietl wykres")
    print("2. Wykonaj odwrotna transformate Fouriera i wyswietl wykres")
    print("3. Proguj sygnal")
    print("4. Zmien plik wejsciowy")
    print("5. Dowolny sinus")
    print("6. Dowolny kosinus")
    print("7. Dowolny sygnal kwadratowy")
    print("8. Dowolny sygnal piloksztaltny")
    print("9. Wyjdz")

def main():
    audio_directory = 'sygnaly' 
    audio_files = list_audio_files(audio_directory)
    
    if not audio_files:
        print("Brak plikow audio w katalogu.")
        return
    
    print("Dostepne pliki audio:")
    for i, file in enumerate(audio_files, 1):
        print(f"{i}. {file}")
    
    file_index = int(input("Wybierz plik audio (numer): ")) - 1
    if file_index < 0 or file_index >= len(audio_files):
        print("Nieprawidłowy wybór.")
        return
    
    selected_file = os.path.join(audio_directory, audio_files[file_index])
    print(f"Wybrano plik: {selected_file}")
    
    while True:
        clear_console()
        display_menu()
        choice = input("Wybierz opcję (numer): ")
        
        if choice == '0':
            print("Wyświetlanie wykresu sygnału audio...")
            m.grafuj(selected_file) 
        elif choice == '1':
            print("Wykonywanie transformaty Fouriera i wyświetlanie wykresu...")
            m.transformata(selected_file) 
        elif choice == '2':
            m.odwrotna(selected_file)
        elif choice == '3':
            prog = float(input("Podaj wartosc progu; "))
            m.progowanie(prog, selected_file)
        elif choice == '4':
            audio_files = list_audio_files(audio_directory)
    
            if not audio_files:
                print("Brak plikow audio w katalogu.")
                return
            print("Dostepne pliki audio:")
            for i, file in enumerate(audio_files, 1):
                print(f"{i}. {file}")
            file_index = int(input("Wybierz plik audio (numer): ")) - 1
            if file_index < 0 or file_index >= len(audio_files):
                print("Nieprawidłowy wybór.")
                return
            selected_file = os.path.join(audio_directory, audio_files[file_index])
            print(f"Wybrano plik: {selected_file}")
        elif choice =='5':
            fprob = 1000
            T = 0.1 
            t = np.linspace(0, T, int(T * fprob), endpoint=False)
            freq = float(input(("Podaj czestotliwosc: ")))
            A = float(input(("Podaj amplitude: ")))
            x = A * np.sin(2 * np.pi * freq * t)
            m.plot(t, x)
        elif choice =='6':
            fprob = 1000
            T = 0.1 
            t = np.linspace(0, T, int(T * fprob), endpoint=False)
            freq = float(input(("Podaj czestotliwosc: ")))
            A = float(input(("Podaj amplitude: ")))
            x = A * np.cos(2 * np.pi * freq * t)
            m.plot(t, x)
        elif choice =='7':
            fprob = 1000
            T = 0.1 
            t = np.linspace(0, T, int(T * fprob), endpoint=False)
            freq = float(input(("Podaj czestotliwosc: ")))
            A = float(input(("Podaj amplitude: ")))
            x = A * signal.square(2 * np.pi * freq * t)
            m.plot(t, x)
        elif choice =='8':
            fprob = 1000
            T = 0.1 
            t = np.linspace(0, T, int(T * fprob), endpoint=False)
            freq = float(input(("Podaj czestotliwosc: ")))
            A = float(input(("Podaj amplitude: ")))
            x = A * signal.sawtooth(2 * np.pi * freq * t)
            m.plot(t, x)      
        elif choice == '9':
            print("Wyjście z programu.")
            break
        else:
            print("Nieprawidłowy wybór, spróbuj ponownie.")

if __name__ == "__main__":
    main()
