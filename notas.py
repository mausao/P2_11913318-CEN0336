#!/usr/bin/env python3
total = 0
contador_notas = 0

while contador_notas < 10:
    nota = input("Digite uma nota: ")
    try: 
        
    except ValueError as e:
        print(f"Erro: {e}" + "\nAs notas devem ser digitos")
    total += float(nota)
    contador_notas += 1

print(f"A média das notas é: ", total/10)
