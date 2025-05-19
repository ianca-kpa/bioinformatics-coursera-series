def LeaderboardCyclopeptideSequencing(Spectrum, N):
    Leaderboard = {""}  # Conjunto contendo apenas o peptídeo vazio
    LeaderPeptides = []  # Lista para armazenar os peptídeos líderes de máximo score
    MaxScore = 0  # Maior score encontrado // MODIFICAÇÃO PARA ESSE EXERCÍCIO
    
    while Leaderboard:
        print(f"Len leaderboard: {len(Leaderboard)}")
        Leaderboard = Expand(Leaderboard)  # Expande os peptídeos no Leaderboard
        for Peptide in list(Leaderboard):  # Itera sobre os peptídeos
            if Mass(Peptide) == ParentMass(Spectrum):
                score = CyclicScore(Peptide, Spectrum, Alphabet, AminoAcidMass)
                if score > MaxScore:
                    MaxScore = score
                    LeaderPeptides = [Peptide]  # Atualiza os peptídeos líderes
                elif score == MaxScore:
                    LeaderPeptides.append(Peptide)
            elif Mass(Peptide) > ParentMass(Spectrum):
                Leaderboard.remove(Peptide)  # Remove peptídeos com massa maior
        Leaderboard = Trim(Leaderboard, Spectrum, N, Alphabet, AminoAcidMass)

    return sorted(LeaderPeptides)  # Garantir saída ordenada

def format_peptide(peptide):
    masses = [str(AminoAcidMass[aa]) for aa in peptide]
    return "-".join(masses)

def Mass(Peptide):
    return sum(AminoAcidMass[aa] for aa in Peptide)

def ParentMass(Spectrum):
    return max(Spectrum)

def Expand(Peptides):
    expanded = []
    for peptide in Peptides:
        for aa in AminoAcidMass.keys():
            expanded.append(peptide + aa)
    return expanded

def CyclicSpectrum(Peptide, Alphabet, AminoAcidMass): # Comentários do código em Bioinfo_II_3
    PrefixMass = [0] * (len(Peptide) + 1)

    for i in range(1, len(Peptide) + 1):
        for s in Alphabet:
            if s == Peptide[i - 1]:
                PrefixMass[i] = PrefixMass[i - 1] + AminoAcidMass[s]

    peptideMass = PrefixMass[len(Peptide)]
    CyclicSpectrum = [0]

    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            CyclicSpectrum.append(PrefixMass[j] - PrefixMass[i])
            if i > 0 and j < len(Peptide):
                CyclicSpectrum.append(
                    peptideMass - (PrefixMass[j] - PrefixMass[i]))

    return sorted(CyclicSpectrum)

def CyclicScore(Peptide, Spectrum, Alphabet, AminoAcidMass):
    Spectrum_Teor = CyclicSpectrum(Peptide, Alphabet, AminoAcidMass)
    score = 0
    i, j = 0, 0
    while i < len(Spectrum_Teor) and j < len(Spectrum):
        if Spectrum_Teor[i] == Spectrum[j]:
            score += 1
            i += 1
            j += 1
        elif Spectrum_Teor[i] < Spectrum[j]:
            i += 1
        else:
            j += 1
    return score

def Trim(Leaderboard, Spectrum, N, Alphabet, AminoAcidMass):
    LinearScores = []
    for peptide in Leaderboard:
        LinearScores.append(LinearScore(peptide, Spectrum, Alphabet, AminoAcidMass))
    sorted_indices = sorted(range(len(LinearScores)), key=lambda x: -LinearScores[x])
    Leaderboard = [Leaderboard[i] for i in sorted_indices]
    LinearScores = [LinearScores[i] for i in sorted_indices]
    if len(Leaderboard) > N:
        cutoff = LinearScores[N - 1]
        Leaderboard = [Leaderboard[i] for i in range(len(Leaderboard)) if LinearScores[i] >= cutoff]
    return Leaderboard

def LinearScore(Peptide, Spectrum, Alphabet, AminoAcidMass):
    Spectrum_Teor = LinearSpectrum(Peptide, Alphabet, AminoAcidMass)
    score = 0
    i, j = 0, 0
    while i < len(Spectrum_Teor) and j < len(Spectrum):
        if Spectrum_Teor[i] == Spectrum[j]:
            score += 1
            i += 1
            j += 1
        elif Spectrum_Teor[i] < Spectrum[j]:
            i += 1
        else:
            j += 1
    return score

def extended_mass_table():
    extended_list = {}
    for i in range(57,201):
        extended_list[chr(i)] = int(i)
    return extended_list

def LinearSpectrum(Peptide, Alphabet, AminoAcidMass):
    PrefixMass = [0] * (len(Peptide) + 1)
    for i in range(1, len(Peptide) + 1):
        for s in Alphabet:
            if s == Peptide[i - 1]: 
                PrefixMass[i] = PrefixMass[i - 1] + AminoAcidMass[s]

    LinearSpectrum = [0]

    for i in range(len(Peptide)):
        for j in range(i + 1, len(Peptide) + 1):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])

    return sorted(LinearSpectrum)

with open("dataset_30245_2.txt", "r") as file:
    lines = file.readlines()
    N = int(lines[0].strip()) 
    Spectrum = list(map(int, lines[1].strip().split())) 

print("aqui")
AminoAcidMass = extended_mass_table()
print("aqui2")
Alphabet = list(AminoAcidMass.keys())
print("aqui3 " + str(N))
result_peptides = LeaderboardCyclopeptideSequencing(Spectrum, N)
print("aqui4")
formatted_result = " ".join(format_peptide(peptide) for peptide in result_peptides)
print(formatted_result)