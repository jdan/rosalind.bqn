AminoOfCodon ⇐ {
  "UUU": 'F'; "UUC": 'F'; "UUA": 'L'; "UUG": 'L';
  "UCU": 'S'; "UCC": 'S'; "UCA": 'S'; "UCG": 'S';
  "UAU": 'Y'; "UAC": 'Y'; "UAA": '*'; "UAG": '*';
  "UGU": 'C'; "UGC": 'C'; "UGA": '*'; "UGG": 'W';
  "CUU": 'L'; "CUC": 'L'; "CUA": 'L'; "CUG": 'L';
  "CCU": 'P'; "CCC": 'P'; "CCA": 'P'; "CCG": 'P';
  "CAU": 'H'; "CAC": 'H'; "CAA": 'Q'; "CAG": 'Q';
  "CGU": 'R'; "CGC": 'R'; "CGA": 'R'; "CGG": 'R';
  "AUU": 'I'; "AUC": 'I'; "AUA": 'I'; "AUG": 'M';
  "ACU": 'T'; "ACC": 'T'; "ACA": 'T'; "ACG": 'T';
  "AAU": 'N'; "AAC": 'N'; "AAA": 'K'; "AAG": 'K';
  "AGU": 'S'; "AGC": 'S'; "AGA": 'R'; "AGG": 'R';
  "GUU": 'V'; "GUC": 'V'; "GUA": 'V'; "GUG": 'V';
  "GCU": 'A'; "GCC": 'A'; "GCA": 'A'; "GCG": 'A';
  "GAU": 'D'; "GAC": 'D'; "GAA": 'E'; "GAG": 'E';
  "GGU": 'G'; "GGC": 'G'; "GGA": 'G'; "GGG": 'G';
  𝕩 : "NotFound"
}

⟨Split⟩ ← •Import "strings.bqn"
# ┌  Fasta ">Rosalind_0498
# │  AAATAAA
# │  >Rosalind_2391
# └  AAATTTT"
# ⟨ ⟨ "Rosalind_0498" "AAATAAA" ⟩ ⟨ "Rosalind_2391" "AAATTT" ⟩ ⟩
Fasta ⇐ {
  Cleanup ← {
    # split by newline
    splitNewLines ← (@+10) Split ¯1↓ 𝕩
    # rejoin as a pair of ⟨LABEL, STRING⟩
    ⟨ (⊑splitNewLines), ∾´1↓splitNewLines ⟩
  }
  Cleanup¨ 1↓ '>' Split 𝕩
}
