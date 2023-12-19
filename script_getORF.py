#!/usr/bin/env python3
import re
import sys
sequencia = {}
codons_frames = {}
protein_frames = {}
translation_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}

#Verificação se o programa está sendo utilizado da forma correta
try:
    if len(sys.argv) != 2:
        print("O uso do programa se dá por: ./script_getORF.py nome_do_arquivo.fasta")
except ValueError as e:
    print(f"Erro: {e}")
    sys.exit(1)

arquivo = sys.argv[1]

#Abertura do arquivo .fasta para manipulação das sequências (agrupamento das sequências em relação ao seu identificador
with open (arquivo) as fasta:
    for linha in fasta:
        linha = linha.rstrip()
        if linha.startswith(">"):
            linhas_juntas = ""
            geneId = re.search(r"^>([\w]{8,10})", linha).group(1)
        else:
            linhas_juntas += linha.strip().upper()
        sequencia[geneId] = "".join(linhas_juntas)

#Manipulação da sequência para obtenção da sequência reversa (calculo dos frames negativos) e separação dos códons por frame
for geneId in sequencia.keys():
    reverseSeq = sequencia[geneId][::-1]
    reverseSeq.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c")
    reverseSeq = reverseSeq.upper()
    codons_frames[geneId] = {}
    for frame in range(3):
        codons = re.findall(r".{3}", sequencia[geneId][frame:])
        frameId = "frame_+" + str(frame+1)
        codons_frames[geneId][frameId] = codons
        codonsReverse = re.findall(r".{3}", reverseSeq[frame:])
        frameId_reverse = "frame_-" + str(frame+1)
        codons_frames[geneId][frameId_reverse] = codonsReverse

#Manipulação dos códons para obtenção da sequência proteica da sequência
for geneId in codons_frames.keys():
    protein_frames[geneId] = {}
    for frame in codons_frames[geneId]:
        protein = ""
        for codon in codons_frames[geneId][frame]:
            if codon in translation_table:
                protein += translation_table[codon]
            else:
                print("O codon não está presente no dicionário de tradução")
                exit()
        protein_frames[geneId][frame] = protein

#Verificação do maior frame aberto de leitura (ORF) e sua respectiva proteina codificada para obtenção da resposta esperada
with open ("ORF.fna", "w") as outputFile:
    with open ("ORF.faa", "w") as outputFile2:
        for geneId in protein_frames.keys():
            longestProtein = ""
            longestFrame = ""
            for frame in protein_frames[geneId]:
                proteinas = re.finditer(r"(M[A-Z]+?\*)", protein_frames[geneId][frame])
            for i in proteinas:
                if len(i.group(1)) > len(longestProtein):
                    longestProtein = i.group(1)
                    longestFrame = frame
                    posStart = i.start(1)
                    posEnd = i.end(1)
                    longestCodon = codons_frames[geneId][frame][posStart:posEnd]
            longestCodon = " ".join(longestCodon)
            resposta = geneId + str(longestFrame) + "_" + str(posStart) + "_" + str(posEnd) + "\n" + longestCodon + "\n"
            outputFile.write(resposta)
            resp = geneId + str(longestFrame) + "_" + str(posStart) + "_" + str(posEnd) + "\n" + longestProtein + "\n"
            outputFile2.write(resp)
