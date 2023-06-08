// Prop
constant symbol Prop : TYPE;
symbol Prf : Prop → TYPE;
constant symbol ⊤ : Prop;
rule Prf ⊤ ↪ Π r, Prf r → Prf r;
constant symbol ⊥ : Prop;
rule Prf ⊥ ↪ Π r, Prf r;
constant symbol ¬ : Prop → Prop; notation ¬ prefix 40;
rule Prf (¬ $p) ↪ Prf $p → Π r, Prf r;
constant symbol ∧ : Prop → Prop → Prop; notation ∧ infix right 30;
rule Prf ($p ∧ $q) ↪ Π r, (Prf $p → Prf $q → Prf r) → Prf r;
constant symbol ∨ : Prop → Prop → Prop; notation ∨ infix right 20;
rule Prf ($p ∨ $q) ↪ Π r, (Prf $p → Prf r) → (Prf $q → Prf r) → Prf r;
constant symbol ⇒ : Prop → Prop → Prop; notation ⇒ infix right 10;
rule Prf ($x ⇒ $y) ↪ Prf $x → Prf $y;

// Set
constant symbol Set : TYPE;
injective symbol El : Set → TYPE;
constant symbol ι : Set;

// Quant
constant symbol ∀ [a : Set] : (El a → Prop) → Prop; notation ∀ quantifier;
rule Prf (∀ $p) ↪ Π x, Prf ($p x);
constant symbol ∃ [a : Set] : (El a → Prop) → Prop; notation ∃ quantifier;
rule Prf (∃ $p) ↪ Π r, (Π x, Prf ($p x) → Prf r) → Prf r;

// Classic
symbol Prfᶜ p ≔ Prf (¬ ¬ p);
symbol ∧ᶜ p q ≔ ¬ ¬ p ∧ ¬ ¬ q; notation ∧ᶜ infix right 30;
symbol ∨ᶜ p q ≔ ¬ ¬ p ∨ ¬ ¬ q; notation ∨ᶜ infix right 20;
symbol ⇒ᶜ p q ≔ ¬ ¬ p ⇒ ¬ ¬ q; notation ⇒ᶜ infix right 10;
symbol ∀ᶜ [a] p ≔ `∀ x : El a, ¬ ¬ (p x); notation ∀ᶜ quantifier;
symbol ∃ᶜ [a] p ≔ `∃ x : El a, ¬ ¬ (p x); notation ∃ᶜ quantifier;
