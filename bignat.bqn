# Arithmetic on large natural numbers.

base ← ⊑ •args ∾ 2⋆24  # Less than half 2⋆53

_while_ ← {𝔽⍟𝔾∘𝔽_𝕣_𝔾∘𝔽⍟𝔾𝕩}

_extend ← { 𝕨 𝔽○((𝕨⌈○≠𝕩)⊸↑) 𝕩 }

Trim ← { 𝕩 ↑˜ -⟜1 _while_ (≥◶⟨0,0=-⊑𝕩˙⟩⟜1) ≠𝕩 }

Carry ← {
  c ← 0
  r ← {𝕩+↩c⋄c↩⌊𝕩÷base⋄base|𝕩}¨ 𝕩
  r‿c
}

Add ⇐ {
  0⊸<◶⊣‿∾´ Carry 𝕨 +_extend 𝕩
}

Sub ⇐ {
  ! 𝕨 ≥○≠ 𝕩
  s‿c ← Carry 𝕨 -_extend 𝕩
  ! 0 = c
  Trim s
}

Mul ⇐ {
  0⊸<◶⟨Trim⊣,∾⟩´ Carry (≠⊸↑˜+⊢)⟜(0⊸∾)˝ 𝕨 ×⌜ 𝕩
}

_compare ⇐ {
  Full ← 0⊸≤◶⟨𝔽˜,⊑⟜𝕨𝔽⊑⟜𝕩⟩ · -⟜1 _while_ (0⊸≤◶⟨0,⊑⟜𝕨=⊑⟜𝕩⟩) -⟜1
  𝕨 =◶𝔽‿Full○≠ 𝕩
}

Expose ⇐ {
  r‿c
}

Neg ⇐ (!0=≠)⊸⊢
Signum ⇐ ×∘≠
Abs ⇐ Conj ⇐ ⊢
