# UBR3 Enrichment Logo Analysis Summary

## Overview
**Analysis Type:** Enrichment-based logo plots (Hits vs Full Screen)  
**Positions Analyzed:** 2-24 (excluding universal Met at position 1)  
**Method:** Height represents frequency × enrichment ratio

---

## 🎯 Key Results

### **Consensus Sequence (Positions 2-24):**
```
M - P D L S V L L L S L S G T G G G S S G G T R L
    │ │ │ │ │ │ │ │ │ │ │ └─────────┘ │ │ │ │ │
    2 3 4 5 6 7 8 9 10... └─G-rich──┘ ... 24
```

---

## 📊 Top 10 Most Discriminating Positions

| Rank | Position | Discrimination Score | Top Enriched AA | % in Hits | Enrichment |
|------|----------|---------------------|-----------------|-----------|------------|
| **1** | **2** | **0.863** | **P (Proline)** | **19.8%** | **3.78x** ⭐ |
| 2 | 12 | 0.714 | M (Met) | 4.4% | 2.31x |
| **3** | **3** | **0.687** | **D (Asp)** | **13.2%** | **2.86x** |
| **4** | **15** | **0.679** | **G (Gly)** | **23.1%** | **2.92x** ⭐⭐ |
| 5 | 6 | 0.678 | V (Val) | 13.2% | 2.22x |
| 6 | 13 | 0.672 | R (Arg) | 14.3% | 2.38x |
| 7 | 22 | 0.629 | I (Ile) | 7.7% | 2.28x |
| 8 | 18 | 0.622 | T (Thr) | 8.8% | 1.86x |
| 9 | 4 | 0.611 | Y (Tyr) | 7.7% | **4.27x** ⭐⭐⭐ |
| 10 | 5 | 0.604 | I (Ile) | 6.6% | 2.04x |

**Note:** Position 4 has the **highest enrichment (4.27x for Tyrosine)** but lower discrimination score due to lower overall frequency.

---

## 🔥 Highest Enrichment Values by Position

### Super-Enriched (>4x):
- **Position 4 - Tyrosine (Y):** 4.27x enrichment (7.7% in hits vs 1.7% in screen)

### Highly Enriched (3-4x):
- **Position 2 - Proline (P):** 3.78x enrichment (19.8% in hits vs 5.2% in screen)
- **Position 23 - Methionine (M):** 3.27x enrichment (5.5% in hits vs 1.6% in screen)

### Moderately Enriched (2-3x):
- **Position 15 - Glycine (G):** 2.92x enrichment (23.1% in hits vs 7.8% in screen)
- **Position 17 - Histidine (H):** 2.92x enrichment (6.6% in hits vs 2.2% in screen)
- **Position 3 - Aspartate (D):** 2.86x enrichment (13.2% in hits vs 4.5% in screen)
- **Position 16 - Methionine (M):** 2.42x enrichment (4.4% in hits vs 1.8% in screen)
- **Position 13 - Arginine (R):** 2.38x enrichment (14.3% in hits vs 5.9% in screen)
- **Position 9 - Histidine (H):** 2.38x enrichment (4.4% in hits vs 1.8% in screen)

---

## 🧬 Functional Regions Identified

### 1. **N-Terminal Recognition Motif (Positions 2-5)**
```
Position 2: P (19.8%, 3.78x) - Turn inducer
Position 3: D (13.2%, 2.86x) - Acidic/charged
Position 4: Y (7.7%, 4.27x) - Aromatic (highest enrichment!)
Position 5: I (6.6%, 2.04x) - Hydrophobic
```

**Interpretation:**
- **Proline at position 2** creates a structural kink
- **Aspartate at position 3** provides negative charge for potential salt bridge
- **Tyrosine at position 4** is the MOST enriched amino acid - may be critical recognition point
- **Hydrophobic at position 5** begins the hydrophobic core

### 2. **Hydrophobic Core (Positions 6-11)**
```
Position 6: V (13.2%, 2.22x)
Position 7: R (15.4%, 2.15x) - Also charged
Position 8: R (15.4%, 2.26x) - Charged
Position 9: H (4.4%, 2.38x) - Charged
```

**Interpretation:**
- Mixed hydrophobic and charged residues
- Arginine/Histidine enrichment suggests electrostatic interactions important
- NOT purely hydrophobic like in full screen

### 3. **Flexible Linker Region (Positions 12-21)** ⭐ CRITICAL
```
Glycine enrichment across this region:
- Position 13: G (14.3%, 1.90x)
- Position 15: G (23.1%, 2.92x) ← HIGHEST FREQUENCY & HIGHLY ENRICHED
- Position 16: G (12.1%, 1.59x)
- Position 17: G (12.1%, 1.53x)
- Position 18: G (13.2%, 1.71x)
- Position 20: G (14.3%, 1.86x)
- Position 21: G (13.2%, 1.61x)
```

**Interpretation:**
- **Position 15 is the hallmark of UBR3 hits** (23.1% glycine, 2.92x enriched)
- Glycine-rich region provides conformational flexibility
- Enables induced-fit binding to UBR3
- This is what distinguishes hits from background library

### 4. **C-Terminal Anchor (Positions 22-24)**
```
Position 22: T/I/N (all ~2x enriched)
Position 23: M (5.5%, 3.27x) - Highly enriched
Position 24: Y/M/Q (all ~1.9x enriched)
```

**Interpretation:**
- Additional contacts with UBR3
- Methionine at position 23 is particularly enriched

---

## 🔬 Detailed Position Analysis

### **Position 2: PROLINE (Most Discriminating)**
- **Frequency:** 19.8% (hits) vs 5.2% (screen)
- **Enrichment:** 3.78x (log2: 1.92)
- **Role:** Turn/kink inducer - constrains backbone conformation
- **Alternative:** Glycine (16.5%, 2.17x) - also flexibility

### **Position 3: ASPARTATE**
- **Frequency:** 13.2% (hits) vs 4.5% (screen)
- **Enrichment:** 2.86x (log2: 1.52)
- **Role:** Acidic residue - potential salt bridge with UBR3
- **Alternative:** Glutamate (13.2%, 1.97x) - also acidic

### **Position 4: TYROSINE** ⭐ HIGHEST ENRICHMENT
- **Frequency:** 7.7% (hits) vs 1.7% (screen)
- **Enrichment:** **4.27x** (log2: 2.09) ← BEST ENRICHMENT
- **Role:** Aromatic ring - may fit into hydrophobic pocket with H-bonding capability
- **Significance:** This is the **most selectively enriched** amino acid

### **Position 15: GLYCINE** ⭐ HIGHEST FREQUENCY
- **Frequency:** **23.1%** (hits) vs 7.8% (screen) ← MOST COMMON IN HITS
- **Enrichment:** 2.92x (log2: 1.54)
- **Role:** Flexibility - allows conformational adaptation
- **Significance:** **Hallmark position** for UBR3 binding

---

## 💡 Biological Insights

### **Proposed UBR3 Recognition Mechanism:**

1. **Initial Recognition (Positions 2-5):**
   - Proline constrains backbone
   - Tyrosine provides key contact (highest specificity)
   - Acidic residue (D/E) forms salt bridge

2. **Binding Interface (Positions 6-11):**
   - Mix of hydrophobic and charged residues
   - Multiple arginine/histidine residues for electrostatic interactions
   - NOT purely hydrophobic binding

3. **Induced Fit (Positions 12-21):**
   - Glycine-rich flexible region
   - Allows substrate to adapt to UBR3 binding groove
   - **Critical discriminating feature**

4. **Anchor Points (Positions 22-24):**
   - Additional contacts for stability
   - Methionine at position 23 particularly important

### **Why This Motif is Selected:**

**Enriched Features:**
- ✅ Structural constraint (Proline at 2)
- ✅ Aromatic contact (Tyrosine at 4) - HIGHEST SPECIFICITY
- ✅ Electrostatic interactions (D/E at 3, R at 7/8/13, H at 9/17)
- ✅ Flexibility (Glycine-rich region 13-21)
- ✅ Specific anchor points (M at 23)

**Avoided Features:**
- ❌ Pure hydrophobic poly-leucine (enriched in background)
- ❌ Alanine (depleted 0.76x)
- ❌ Large aliphatic (Leucine depleted 0.81x)

---

## 📈 Files Generated

### **Primary Results:**
1. **`enrichment_logo_comprehensive.png`** ⭐⭐⭐
   - 5-panel comparison showing:
     - Hits sequence logo
     - Screen sequence logo
     - **Enrichment-weighted logo (KEY PLOT)**
     - Log2 enrichment heatmap
     - Position discrimination scores

2. **`enrichment_logo_simplified.png`** ⭐⭐
   - Clean enrichment logo only
   - Best for presentations

3. **`position_enrichment_plots.png`** ⭐
   - Top 8 discriminating positions in detail
   - Bar plots comparing hits vs screen

### **Data Files:**
- `enrichment_matrix_pos2+.csv` - Raw enrichment ratios
- `log2_enrichment_pos2+.csv` - Log2 transformed enrichment
- `enrichment_weighted_pos2+.csv` - Frequency × enrichment matrix
- `enrichment_summary.txt` - Complete text summary

---

## 🎯 Key Takeaways

### **Top 3 Critical Positions:**
1. **Position 2 (Proline)** - Most discriminating overall (3.78x)
2. **Position 4 (Tyrosine)** - Highest enrichment (4.27x) ⭐
3. **Position 15 (Glycine)** - Hallmark of hits (23.1% frequency, 2.92x)

### **Recognition Signature:**
```
M-[P/G]-[D/E]-[Y]-[hydrophobic]-[charged]-[FLEXIBLE GLYCINE-RICH]-[anchor]
│   │     │    │                            └─positions 12-21─┘
1   2     3    4 (HIGHEST ENRICHMENT)                ↑
                                            MOST DISTINCTIVE FEATURE
```

### **Novel Findings:**
1. **Tyrosine at position 4** shows the highest fold-enrichment (4.27x) - previously unappreciated
2. **Glycine at position 15** is the most abundant discriminating feature (23.1%)
3. **Mixed hydrophobic/charged character** distinguishes from typical hydrophobic degrons
4. **Flexibility is critical** - glycine enrichment across positions 13-21

---

## 🔬 Recommended Experiments

### High Priority:
1. **Mutagenesis of position 4 (Y):**
   - Y4A (remove aromatic)
   - Y4F (keep aromatic, remove OH)
   - Y4W (larger aromatic)
   - **Expected:** Loss of UBR3 binding with Y4A

2. **Mutagenesis of position 2 (P):**
   - P2A (remove kink)
   - P2G (flexible)
   - **Expected:** Reduced but not abolished binding

3. **Glycine region modification:**
   - Replace G15 with A
   - Replace multiple Gs (13-21) with A
   - **Expected:** Major loss of binding due to reduced flexibility

### Structural Studies:
- Co-crystallize UBR3 + peptide with consensus motif
- Focus on Y4 interaction
- Measure flexibility of G-rich region by NMR

---

**Analysis Date:** January 17, 2026  
**Method:** Enrichment-based logo plots (Hits/Screen ratios)  
**Sequences:** 91 hits vs 16,514 full screen  
**Positions:** 2-24 (excluding universal M1)
