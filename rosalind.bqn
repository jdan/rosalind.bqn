⟨ Fasta, AminoOfCodon ⟩ ← •Import "utils.bqn"

DNA ← +´=⟜"ACGT"¨
DNA "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

DNA2 ← +˝"ACGT"⊸(=⌜˜) # can also top-left = and and ˝ insert
DNA2 "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

RNA ← {'T':'U';𝕩}¨
RNA "GATGGAACTTGACTACGTAAATT"

REVC ← ⌽{'A':'T';'T':'A';'C':'G';'G':'C'}¨
REVC "AAAACCCGGT"

FIB ← { ⟨n,k⟩:
  # a, b -> b, b+k*a
  Next ← { ⟨ 1⊑𝕩, 1⊑𝕩 + k×⊑𝕩 ⟩ }
  ⊑ Next⍟n 0‿1
}
FIB 5‿3

GC ← {𝕊 input:
  # TODO: Parses FASTA format, used in other problems
  parsed ← Fasta input

  # Compute the GC content of a DNA string
  GCpct ← {𝕊 dna:
    gc ← +´ { (𝕩 = 'G') ∨ (𝕩 = 'C') }¨ dna
    100×gc÷≠dna
  }

  # Return the label and GC content of the DNA string
  # with the highest GC content
  labels ← ⊑¨ parsed
  gcs ← GCpct∘(1⊸⊑)¨ parsed
  maxIndex ← ⊑ ⍒ gcs

  maxIndex⊑¨ ⟨ labels, gcs ⟩
}

GC ">Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT
"

HAMM ← {⟨a,b⟩: +´a≠b}
HAMM "GAGCCTACTAACGGGAT"‿"CATCGTAATGACGGCCT"

IPRB ← {⟨k,m,n⟩:
  +´ ⟨
    # P(k𝕩) = 1
    k÷(k+m+n),
    # P(mk) = 1
    (m÷(k+m+n))×k÷(k+m+n-1),
    # P(mm) = 3/4
    (m÷(k+m+n))×((m-1)÷(k+m+n-1))×3÷4,
    # P(mn) = 1/2
    (m÷(k+m+n))×(n÷(k+m+n-1))×1÷2,
    # P(nk) = 1
    (n÷(k+m+n))×k÷(k+m+n-1),
    # P(nm) = 1/2
    (n÷(k+m+n))×(m÷(k+m+n-1))×1÷2
  ⟩
}
IPRB 2‿2‿2

SUBS ← { ⟨needle,haystack⟩:
  len ← ≠ haystack
  matches ← needle ⍷ haystack
  # indices of matches
  1 + (len ↑ matches) / ↕len
}
SUBS ⟨"ATAT", "GATATATGCATATACTT"⟩

# TODO: Need double-precision numbers, bignat
# likely not the right choice.
FIBD ← {⟨n,m⟩:
  ⟨Add⟩ ← •Import "bignat.bqn"
  zero ← 0 Add 0
  one ← 0 Add 1

  Next ← {𝕊 rabbits:
    newRabbits ← Add´ 1↓rabbits
    # drop the last group
    ⟨newRabbits⟩ ∾ (¯1↓rabbits)
  }
  # start with a single pair of 0-month
  # old rabbits, and complete the state
  # with m-1 empty age groups
  init ← ⟨⟨one⟩⟩ ∾ (m-1) / ⋈⋈zero
  ⊑ Add´ Next⍟(n-1) init
}
FIBD 52‿20

# P(AA × AA) = 1
# P(AA × Aa) = 1
# P(AA × aa) = 1
# P(Aa × Aa) = 3/4
# P(Aa × aa) = 1/2
# P(aa × aa) = 0
IEV ← {⟨a,b,c,d,e,f⟩:(2×(a+b+c))+(1.5×d)+e}
IEV 18916‿18455‿17826‿16977‿19645‿19837  # 155504.5

LIA ← {⟨k,n⟩:
  # The key to this problem is realizing  that no matter what genotype the
  # organism that mates with someone of genotype Aa Bb has, the probability that
  # the offspring is Aa Bb is always 0.25
  p ← 1÷4
  Fact ← {×´(1+↕𝕩)}
  Choose ← {(Fact 𝕨) ÷ ((Fact 𝕩) × (Fact 𝕨 - 𝕩))}
  # P(r) = \binom{2^k}{r} \cdot p^r \cdot (1-p)^{2^k-r}
  PExact ← {((2⋆k) Choose 𝕩) × (p⋆𝕩) × ((1-p)⋆((2⋆k) - 𝕩))}
  Range ← {𝕨 + ↕(1 + 𝕩 - 𝕨)}
  # P(r≥n) = \sum_{i=n}^{2^k} P(i)
  +´ PExact¨ n Range 2⋆k
}
LIA 5‿7

PRTM ← {𝕊 protein:
  MonoisotopicMass ← {
    'A': 71.03711;'C': 103.00919;'D': 115.02694;
    'E': 129.04259;'F': 147.06841;'G': 57.02146;
    'H': 137.05891;'I': 113.08406;'K': 128.09496;
    'L': 113.08406;'M': 131.04049;'N': 114.04293;
    'P': 97.05276;'Q': 128.05858;'R': 156.10111;
    'S': 87.03203;'T': 101.04768;'V': 99.06841;
    'W': 186.07931;'Y': 163.06333
  }
  +´ MonoisotopicMass¨ protein
}
PRTM "SKADYEK"

PROT ← {
  ⟨⟩: ⟨⟩ ;
  # TODO: Stop recursing at "*"
  𝕩: (AminoOfCodon 3↑𝕩) ∾ 𝕊 3↓𝕩
}
PROT "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"

GRPH ← {⟨k, input⟩:
  parsed ← Fasta input
  dnas ← 1⊑¨ parsed
  Matches ← {
    (𝕩 ≡ 𝕨) ? 0 ;
    ((- k)↑ 𝕨) ≡ (k↑𝕩)
  }
  # adjacency table of suffix = prefix
  adj ← Matches⌜˜ dnas
  # turn the table into a list of edges
  MatchesForEntry ← { 𝕊 n:
    others ← (n ⊏ adj) / ⊑¨ parsed
    (⊑n⊑parsed) ⊸⋈¨ others
  }
  PrintLine ← {⟨label, dna⟩:
    nl ← @ + 10
    label ∾ " " ∾ dna ∾ nl
  }
  ∾´ PrintLine¨ ∾´ MatchesForEntry¨ ↕≠parsed
}

GRPH 3‿">Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG
"
