# -*- coding: utf-8 -*-
# backend/mirna_engine.py -- NanoSynapse Engine v3.0 (User-Toggleable Biological Rules)
# v3.0 changes:
#   [1] 5 biological rules are now individually user-toggleable via parameters
#   [2] Rule iii (pos 1-9 max 1 mismatch) — newly implemented (was missing)
#   [3] Rule iv (consecutive 2+ ban) — reset logic bug fixed
#   [4] Rule v  (seed region MFE ≤ threshold) — calculate_region_mfe() added
#   [5] evaluate_strict_rules() signature extended with 5 rule flags
#   [6] find_targets() signature extended with rule flags + seed_mfe_threshold
#   [7] Animal mode cleavage override removed — now fully user-controlled

import random
import RNA
import numpy as np
from Bio.SeqUtils import gc_fraction

# ─── CONSTANTS ────────────────────────────────────────────────────────────────

# Position-based weights for Animal Mode (0-based indexes for miRNA positions 1-22)
ANIMAL_POSITION_WEIGHTS = {
    0: 0.5,
    1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0,  # Core Seed (pos 2-8)
    8: 0.5,
    9: 1.5, 10: 1.5,   # Argonaute RISC cleavage site (pos 10-11)
    11: 1.0, 12: 1.0, 13: 1.0, 14: 1.0, 15: 1.0,
    16: 0.5, 17: 0.5, 18: 0.5, 19: 0.5, 20: 0.5, 21: 0.5
}

WATSON_CRICK = [frozenset({'A', 'U'}), frozenset({'C', 'G'})]
WOBBLE       = frozenset({'G', 'U'})
ALL_PAIRS    = WATSON_CRICK + [WOBBLE]

DEFAULT_MAX_PENALTY_PLANT   = 4.0   # psRNATarget standard (was 4.5)
DEFAULT_MAX_PENALTY_ANIMAL  = 4.0   # Miranda-equivalent (was 5.0, hardcoded)
DEFAULT_MFE_THRESHOLD_ANIMAL = -17.0  # Relaxed to reduce false negatives (was -20.0)
DEFAULT_MFE_THRESHOLD_PLANT  = -18.0  # Tightened to reduce false positives (was -15.0)

# ─── 1. ALIGNMENT PARSING ─────────────────────────────────────────────────────

def parse_alignment(mirna_seq, target_seq, structure, duplex_i, duplex_j):
    """
    Correctly parse ViennaRNA dot-bracket structure and extract exact alignment columns.
    Maps anti-parallel strands (miRNA 5'->3' and target 3'->5') including bulges/gaps.
    
    mirna_seq: 5' to 3' miRNA sequence
    target_seq: 5' to 3' full target sequence
    structure: dot-bracket string from duplexfold (e.g. '(((...)))&(((...)))')
    duplex_i: 1-based end index in miRNA
    duplex_j: 1-based start index in target
    """
    if not structure or '&' not in structure:
        return []
        
    s1, s2 = structure.split('&', 1)
    
    # Extract aligned portion of miRNA
    mirna_len = len(s1)
    mirna_start = max(0, duplex_i - mirna_len)
    mirna_aligned = mirna_seq[mirna_start:duplex_i]
    
    # Extract aligned portion of target
    target_len = len(s2)
    target_start = max(0, duplex_j - 1)
    target_end = min(len(target_seq), duplex_j - 1 + target_len)
    target_aligned = target_seq[target_start:target_end]
    
    # Find base pairing indices (representing Watson-Crick and Wobble pairs)
    s1_pairs = [idx for idx, char in enumerate(s1) if char == '(']
    s2_pairs = [idx for idx, char in enumerate(s2) if char == ')']
    
    # Reverse target pairing since the hybrid binds anti-parallel
    s2_pairs_reversed = list(reversed(s2_pairs))
    
    pairing_map = {}
    for k in range(min(len(s1_pairs), len(s2_pairs))):
        pairing_map[s1_pairs[k]] = s2_pairs_reversed[k]
        
    columns = []
    target_placed = [False] * len(target_aligned)
    
    # Loop over miRNA position i from 5' to 3'
    for i in range(len(s1)):
        char1 = s1[i]
        
        if char1 == '(':
            target_idx = pairing_map[i]
            
            # Place any target bulges before this paired base.
            # Antiparallel check: from the right (end) of target_aligned down to target_idx + 1.
            for j in range(len(target_aligned) - 1, target_idx, -1):
                if j < len(target_placed) and not target_placed[j]:
                    columns.append({
                        "m": "-",
                        "t": target_aligned[j],
                        "match": False,
                        "bulge": "target"
                    })
                    target_placed[j] = True
            
            # Place the paired column
            if i < len(mirna_aligned) and target_idx < len(target_aligned):
                m_base = mirna_aligned[i]
                t_base = target_aligned[target_idx]
                pair = frozenset({m_base, t_base})
                is_wc = pair in WATSON_CRICK
                is_wob = pair == WOBBLE
                
                columns.append({
                    "m": m_base,
                    "t": t_base,
                    "match": is_wc or is_wob
                })
                target_placed[target_idx] = True
        else:
            # miRNA bulge (gap in target)
            if i < len(mirna_aligned):
                columns.append({
                    "m": mirna_aligned[i],
                    "t": "-",
                    "match": False,
                    "bulge": "mirna"
                })
                
    # Place any remaining target bulges at the very 5' end of the target
    # (which corresponds to the 3' end of the miRNA)
    for j in range(len(target_aligned) - 1, -1, -1):
        if j < len(target_placed) and not target_placed[j]:
            columns.append({
                "m": "-",
                "t": target_aligned[j],
                "match": False,
                "bulge": "target"
            })
            target_placed[j] = True
            
    return columns

# ─── 2. WEIGHTED PENALTY SCORE ────────────────────────────────────────────────

def calculate_penalty_score(alignment, rule_set="animal"):
    """
    Weighted mismatch and gap penalty score.
    
    For Plant (psRNATarget expectation scoring):
      - Perfect WC match = 0 penalty.
      - Wobble G:U = 0.5.
      - Mismatch = 1.0.
      - Mismatches in seed (positions 2-13) are multiplied by 1.5. Wobbles do NOT get multiplier.
      - Gap opening = 2.0.
      - Gap extension = 0.5.
      
    For Animal:
      - Weighted mismatches/wobbles using ANIMAL_POSITION_WEIGHTS.
      - Gap opening = 2.0 * position weight.
      - Gap extension = 0.5 * position weight.
    """
    score = 0.0
    
    if rule_set == "plant":
        mirna_pos = 0
        in_gap = False
        
        for col in alignment:
            m = col["m"]
            t = col["t"]
            
            if m == "-":
                # Target bulge (gap in miRNA)
                if not in_gap:
                    score += 2.0  # Gap opening
                    in_gap = True
                else:
                    score += 0.5  # Gap extension
                continue
                
            if t == "-":
                # miRNA bulge (gap in target)
                if not in_gap:
                    score += 2.0  # Gap opening
                    in_gap = True
                else:
                    score += 0.5  # Gap extension
                mirna_pos += 1
                continue
                
            in_gap = False
            mirna_pos += 1  # 1-based miRNA position
            
            pair = frozenset({m, t})
            if pair in WATSON_CRICK:
                pass
            elif pair == WOBBLE:
                score += 0.5  # Wobble
            else:
                # Mismatch - check seed multiplier (miRNA pos 2-13)
                is_seed = (2 <= mirna_pos <= 13)
                multiplier = 1.5 if is_seed else 1.0
                score += 1.0 * multiplier
                
    else:
        # Animal mode (using ANIMAL_POSITION_WEIGHTS)
        mirna_pos = 0
        in_gap = False
        
        for col in alignment:
            m = col["m"]
            t = col["t"]
            
            if m == "-":
                w = ANIMAL_POSITION_WEIGHTS.get(mirna_pos, 1.0)
                if not in_gap:
                    score += 2.0 * w
                    in_gap = True
                else:
                    score += 0.5 * w
                continue
                
            if t == "-":
                w = ANIMAL_POSITION_WEIGHTS.get(mirna_pos, 1.0)
                if not in_gap:
                    score += 2.0 * w
                    in_gap = True
                else:
                    score += 0.5 * w
                mirna_pos += 1
                continue
                
            in_gap = False
            w = ANIMAL_POSITION_WEIGHTS.get(mirna_pos, 1.0)
            mirna_pos += 1
            
            pair = frozenset({m, t})
            if pair in WATSON_CRICK:
                pass
            elif pair == WOBBLE:
                score += 0.5 * w
            else:
                score += 1.0 * w
                
    return round(score, 2)

# ─── 3. STRICT VALIDATION RULES ───────────────────────────────────────────────

def evaluate_strict_rules(
    alignment,
    rule_set="animal",
    max_mismatches=4,
    strict_cleavage=True,
    # ── Granüler kural toggleları (v3.0) ──────────────────────────────────────
    rule_max_mismatch=True,   # Kural i:   Global max mismatch kontrolü
    rule_cleavage_site=True,  # Kural ii:  Pos 10-11 cleavage koruması
    rule_seed_mismatch=True,  # Kural iii: Pos 1-9 max 1 mismatch
    rule_consecutive=True,    # Kural iv:  Ardışık 2+ mismatch/gap yasak
):
    """
    Evaluates strict alignment validation rules.
    Each of the 5 biological rules can be individually toggled.
    Returns (passed, reasons, mismatches, gaps, penalty)
    
    Rules:
      i.   rule_max_mismatch  — Total mismatches+gaps <= max_mismatches
      ii.  rule_cleavage_site — No mismatch at pos 10-11 (5' end)
      iii. rule_seed_mismatch — Max 1 mismatch in pos 1-9
      iv.  rule_consecutive   — No more than 2 consecutive mismatches/gaps
      (v handled in find_targets via calculate_region_mfe)
    """
    passed = True
    reasons = []
    mismatches = 0
    gaps = 0
    consecutive = 0
    max_consecutive_seen = 0
    seed_mm = 0
    seed_wobbles = 0
    early_seed_mm = 0   # Kural iii: pos 1-9 mismatch sayacı

    mirna_pos = 0

    for col in alignment:
        m = col["m"]
        t = col["t"]

        # ── GAP sütunu ──────────────────────────────────────────────────────────
        if m == "-" or t == "-":
            gaps += 1
            consecutive += 1
            if consecutive > max_consecutive_seen:
                max_consecutive_seen = consecutive

            if m != "-":  # miRNA gap (bulge in target)
                mirna_pos += 1
                if rule_set == "animal" and (2 <= mirna_pos <= 8):
                    seed_mm += 1
                elif rule_set == "plant" and (2 <= mirna_pos <= 13):
                    seed_mm += 1
                # Kural iii: pos 1-9 gap sayılır
                if 1 <= mirna_pos <= 9:
                    early_seed_mm += 1
            continue

        # ── PAIRED sütunu ───────────────────────────────────────────────────────
        mirna_pos += 1
        pair = frozenset({m, t})
        is_match = pair in ALL_PAIRS
        is_wobble = pair == WOBBLE

        if is_match:
            # Düzgün çift: consecutive sıfırla
            consecutive = 0
            if is_wobble:
                if rule_set == "animal" and (2 <= mirna_pos <= 8):
                    seed_wobbles += 1
                elif rule_set == "plant" and (2 <= mirna_pos <= 13):
                    seed_wobbles += 1
        else:
            # Mismatch
            mismatches += 1
            consecutive += 1
            if consecutive > max_consecutive_seen:
                max_consecutive_seen = consecutive

            # Seed region sayacı
            if rule_set == "animal" and (2 <= mirna_pos <= 8):
                seed_mm += 1
            elif rule_set == "plant" and (2 <= mirna_pos <= 13):
                seed_mm += 1

            # Kural iii: pos 1-9 mismatch sayacı
            if 1 <= mirna_pos <= 9:
                early_seed_mm += 1

            # Kural ii: Cleavage site koruması (pos 10-11)
            if rule_cleavage_site and strict_cleavage and (10 <= mirna_pos <= 11):
                passed = False
                reasons.append("Kural ii: Mismatch at cleavage site pos {} (pos 10-11 protected)".format(mirna_pos))

    # ── Kural iv: Ardışık mismatch/gap kontrolü ──────────────────────────────
    if rule_consecutive and max_consecutive_seen > 2:
        passed = False
        reasons.append(
            "Kural iv: {} consecutive mismatches/gaps found (max 2 allowed)".format(max_consecutive_seen)
        )

    # ── Kural iii: Pos 1-9 max 1 mismatch ───────────────────────────────────
    if rule_seed_mismatch and early_seed_mm > 1:
        passed = False
        reasons.append(
            "Kural iii: {} mismatches/gaps in pos 1-9 (max 1 allowed)".format(early_seed_mm)
        )

    # ── Rule-set spesifik seed kontrolleri ───────────────────────────────────
    if rule_set == "animal":
        if seed_mm > 2:
            passed = False
            reasons.append("Too many seed mismatches/gaps in pos 2-8 ({} found, max 2)".format(seed_mm))
        if seed_wobbles > 2:
            passed = False
            reasons.append("Too many G:U wobbles in animal seed pos 2-8 ({} found, max 2)".format(seed_wobbles))

        # Kural i: Global total mismatch kontrolü
        if rule_max_mismatch:
            total_errors = mismatches + gaps
            if total_errors > max_mismatches:
                passed = False
                reasons.append(
                    "Kural i: Total mismatches & gaps ({}) exceeds limit {}".format(total_errors, max_mismatches)
                )
    else:
        # Plant: seed 2-13
        if seed_mm > 2:
            passed = False
            reasons.append("Too many seed mismatches/gaps in pos 2-13 ({} found, max 2)".format(seed_mm))

        # Plant seed'de ardışık mismatch
        consec_seed_mm = 0
        mirna_pos_temp = 0
        for col in alignment:
            if col["m"] == "-":
                continue
            mirna_pos_temp += 1
            if 2 <= mirna_pos_temp <= 13:
                if col["t"] == "-" or (frozenset({col["m"], col["t"]}) not in ALL_PAIRS):
                    consec_seed_mm += 1
                    if consec_seed_mm > 1:
                        passed = False
                        reasons.append("Consecutive mismatches/gaps in seed region not allowed")
                        break
                else:
                    consec_seed_mm = 0

        # Kural i: Global total mismatch kontrolü (plant)
        if rule_max_mismatch:
            total_errors = mismatches + gaps
            if total_errors > max_mismatches:
                passed = False
                reasons.append(
                    "Kural i: Total mismatches & gaps ({}) exceeds limit {}".format(total_errors, max_mismatches)
                )

    # ── Penalty score değerlendirmesi ────────────────────────────────────────
    penalty = calculate_penalty_score(alignment, rule_set)
    max_penalty = DEFAULT_MAX_PENALTY_PLANT if rule_set == "plant" else DEFAULT_MAX_PENALTY_ANIMAL
    if penalty > max_penalty:
        passed = False
        reasons.append("Penalty score {} exceeds threshold {}".format(penalty, max_penalty))

    return passed, list(set(reasons)), mismatches, gaps, penalty

# ─── 4. SEED REGION MFE (Kural v) ────────────────────────────────────────────

def calculate_region_mfe(mirna_seq, target_seq, bind_start):
    """
    Kural v: miRNA 5' ucundan pos 2-7 ve/veya pos 8-13 için bölge-spesifik MFE.
    Her iki seed penceresi için ayrı duplexfold hesaplar.
    En iyi (en negatif) MFE değerini döndürür.
    
    mirna_seq:  tam miRNA dizisi (5'→3')
    target_seq: hedef sekans
    bind_start: 1-tabanlı bağlanma başlangıcı
    """
    try:
        seed_2_7  = mirna_seq[1:7]    # 0-indexed → pos 2-7
        seed_8_13 = mirna_seq[7:13]   # 0-indexed → pos 8-13

        # Hedef üzerinde bağlanma bölgesi etrafında pencere
        ctx_start = max(0, bind_start - 5)
        ctx_end   = min(len(target_seq), bind_start + 20)
        window    = target_seq[ctx_start:ctx_end]

        if not window or len(window) < 6:
            return 0.0, 0.0

        mfe_2_7  = RNA.duplexfold(seed_2_7,  window).energy if len(seed_2_7)  >= 4 else 0.0
        mfe_8_13 = RNA.duplexfold(seed_8_13, window).energy if len(seed_8_13) >= 4 else 0.0

        return round(mfe_2_7, 2), round(mfe_8_13, 2)
    except Exception:
        return 0.0, 0.0


# ─── 5. TARGET SITE ACCESSIBILITY ─────────────────────────────────────────────

def estimate_accessibility(target_seq, bind_start, mirna_len):
    """
    Fraction of unpaired bases in the binding window.
    Uses RNA.fold on a ±20 nt context window around the binding site.
    """
    try:
        # bind_start is 1-based index in target sequence
        ctx_start = max(0, bind_start - 21)
        ctx_end   = min(len(target_seq), bind_start - 1 + mirna_len + 20)
        region    = target_seq[ctx_start:ctx_end]
        
        structure, _ = RNA.fold(region)
        
        rel_start = bind_start - 1 - ctx_start
        rel_end   = rel_start + mirna_len
        window    = structure[rel_start : min(rel_end, len(structure))]
        
        if not window:
            return 1.0
            
        paired = window.count('(') + window.count(')')
        return round(1.0 - paired / len(window), 3)
    except Exception:
        return 1.0  # Default to accessible on error

# ─── 5. ALL BINDING SITES (SLIDING WINDOW) ────────────────────────────────────

def find_all_binding_sites(mirna_seq, target_seq, min_mfe=-15.0, max_sites=3):
    """
    Slides a window across the target sequence to find thermodynamic binding sites.
    Returns sorted top-N sites with correct coordinates.
    """
    win = len(mirna_seq) + 12
    step = 5
    sites = []
    
    for i in range(0, max(1, len(target_seq) - win + 1), step):
        window = target_seq[i : i + win]
        duplex = RNA.duplexfold(mirna_seq, window)
        
        if duplex.energy >= min_mfe:
            continue
            
        s1, s2 = duplex.structure.split('&')
        exact_start = i + duplex.j - 1
        exact_end = exact_start + len(s2)
        
        merged = False
        for site in sites:
            if abs(exact_start - site['start']) < win // 2:
                if duplex.energy < site['mfe']:
                    site.update(start=exact_start, end=exact_end,
                                mfe=round(duplex.energy, 2),
                                structure=duplex.structure)
                merged = True
                break
                
        if not merged:
            sites.append(dict(start=exact_start, end=exact_end,
                              mfe=round(duplex.energy, 2),
                              structure=duplex.structure))
                              
    return sorted(sites, key=lambda x: x['mfe'])[:max_sites]

# ─── 6. P-VALUE (MONTE CARLO PERMUTATION TEST) ────────────────────────────────

def calculate_pvalue(mirna_seq, target_seq, actual_mfe, n=200):  # Increased from 50 → 200 for statistical reliability
    """
    Shuffle miRNA n times, compute background MFE distribution.
    p-value = fraction of shuffled alignments with MFE <= actual_mfe.
    """
    try:
        seq_list = list(mirna_seq)
        background = []
        target_window = target_seq[:120]  # Use 120 nt window for speed
        
        for _ in range(n):
            random.shuffle(seq_list)
            duplex = RNA.duplexfold(''.join(seq_list), target_window)
            background.append(duplex.energy)
            
        count_lte = sum(1 for m in background if m <= actual_mfe)
        return round(max(float(count_lte) / n, 1.0 / n), 4)
    except Exception:
        return 1.0

# ─── 7. COMPOSITE CONFIDENCE SCORE ────────────────────────────────────────────

def calculate_confidence(mfe, penalty, gc, accessibility, p_value, rule_set="animal"):
    """
    Composite score 0-100 combining all biological features.
    
    For animals: Focuses heavily on thermodynamic energy (MFE).
      Weights: MFE 45% | Penalty 25% | Accessibility 15% | GC 10% | p-value 5%
    For plants: Focuses heavily on structural alignment penalty.
      Weights: MFE 25% | Penalty 45% | Accessibility 15% | GC 10% | p-value 5%
    # v2.2 weights: penalty importance increased, p-value strengthened, GC reduced
    # Animal: MFE 40% | Penalty 30% | Accessibility 15% | GC 5% | p-value 10%
    # Plant:  MFE 25% | Penalty 45% | Accessibility 15% | GC 5% | p-value 10%
    """
    mfe_s    = min(100.0, abs(mfe) * 2.5)           # -40 kcal → 100
    pen_s    = max(0.0,   100.0 - penalty * 20)     #  0 penalty  → 100
    acc_s    = accessibility * 100.0                 #  1.0      → 100
    gc_s     = max(0.0,   100.0 - abs(gc - 50) * 2) # GC=50%   → 100
    pval_s   = max(0.0,   (1.0 - p_value) * 100)    # p=0       → 100
    
    if rule_set == "plant":
        val = mfe_s * 0.25 + pen_s * 0.45 + acc_s * 0.15 + gc_s * 0.05 + pval_s * 0.10
    else:
        val = mfe_s * 0.40 + pen_s * 0.30 + acc_s * 0.15 + gc_s * 0.05 + pval_s * 0.10
        
    return round(min(100.0, val), 1)

# ─── MAIN PREDICTION FUNCTION ─────────────────────────────────────────────────

def find_targets(
    mirna_id, mirna_seq, target_id, target_seq,
    compute_pvalue=False,
    rule_set="animal",
    max_mismatches=4,
    strict_cleavage=True,
    mfe_threshold=None,
    # ── v3.0: Granüler kural toggleları ──────────────────────────────────────
    rule_max_mismatch=True,   # Kural i
    rule_cleavage_site=True,  # Kural ii
    rule_seed_mismatch=True,  # Kural iii
    rule_consecutive=True,    # Kural iv
    rule_seed_mfe=False,      # Kural v (default kapalı, ek hesaplama)
    seed_mfe_threshold=-20.0, # Kural v eşiği
):
    """
    Full miRNA-target prediction pipeline supporting plant and animal rule sets.
    Returns a rich result dictionary consumed by FastAPI and Next.js frontend.
    """
    mirna_len = len(mirna_seq)
    
    # Step 1: MFE via ViennaRNA duplexfold
    duplex = RNA.duplexfold(mirna_seq, target_seq)
    mfe = duplex.energy
    structure = duplex.structure
    
    # Split structure to extract correct coordinates
    s1, s2 = structure.split('&', 1)
    bind_start = duplex.j
    bind_end = duplex.j + len(s2) - 1
    
    # Slice the correct matched target subsequence
    slice_start = max(0, bind_start - 1)
    slice_end = min(len(target_seq), bind_end)
    matched = target_seq[slice_start:slice_end]
    
    # Step 2: Parse exact alignment BasePair[] columns for frontend
    alignment = parse_alignment(mirna_seq, target_seq, structure, duplex.i, duplex.j)
    
    # Step 3: Evaluate strict rules and weighted penalty score along alignment columns
    passed, reasons, mm_count, gaps, penalty = evaluate_strict_rules(
        alignment, rule_set, max_mismatches, strict_cleavage,
        rule_max_mismatch=rule_max_mismatch,
        rule_cleavage_site=rule_cleavage_site,
        rule_seed_mismatch=rule_seed_mismatch,
        rule_consecutive=rule_consecutive,
    )
    
    # Step 4: Check thermodynamic threshold
    threshold_mfe = mfe_threshold if mfe_threshold is not None else (DEFAULT_MFE_THRESHOLD_PLANT if rule_set == "plant" else DEFAULT_MFE_THRESHOLD_ANIMAL)
    if mfe > threshold_mfe:
        passed = False
        reasons.append("MFE ({:.1f} kcal/mol) above threshold {}".format(mfe, threshold_mfe))
        
    # Step 5: Target site accessibility
    accessibility = estimate_accessibility(target_seq, bind_start, mirna_len)
    
    # Step 6: GC content
    gc_pct = round(gc_fraction(mirna_seq) * 100, 2)
    
    # Step 7: p-value (only for manual alignments due to db scan speeds)
    p_value = calculate_pvalue(mirna_seq, target_seq, mfe) if compute_pvalue else 1.0

    # Step 7b: Kural v — Seed region MFE (pos 2-7 ve 8-13) — sadece gerektiğinde
    seed_mfe_2_7  = 0.0
    seed_mfe_8_13 = 0.0
    if rule_seed_mfe:
        seed_mfe_2_7, seed_mfe_8_13 = calculate_region_mfe(mirna_seq, target_seq, bind_start)
        if not (seed_mfe_2_7 <= seed_mfe_threshold or seed_mfe_8_13 <= seed_mfe_threshold):
            passed = False
            reasons.append(
                "Kural v: Seed MFE insufficient — pos2-7: {} kcal/mol, pos8-13: {} kcal/mol (need <= {} in at least one)".format(
                    seed_mfe_2_7, seed_mfe_8_13, seed_mfe_threshold
                )
            )

    # Step 8: All binding sites
    min_site_mfe = DEFAULT_MFE_THRESHOLD_PLANT if rule_set == "plant" else DEFAULT_MFE_THRESHOLD_ANIMAL
    all_sites = find_all_binding_sites(mirna_seq, target_seq, min_site_mfe) if len(target_seq) > mirna_len + 12 else []
    
    # Step 9: Composite confidence score
    confidence = calculate_confidence(mfe, penalty, gc_pct, accessibility, p_value, rule_set)
    
    # Step 10: Similarity percentage (matches out of total columns)
    aligned_len = len(alignment)
    match_count = sum(1 for col in alignment if col["match"])
    similarity = round((float(match_count) / max(1, aligned_len)) * 100, 2)
    
    return {
        "miRNA_ID":           mirna_id,
        "Gene_ID":            target_id,
        "Status":             "PASS" if passed else "FAIL",
        "Confidence_Score":   confidence,
        "MFE_kcal_mol":       round(mfe, 2),
        "Penalty_Score":      penalty,
        "Accessibility":      accessibility,
        "GC_Content_Percent": gc_pct,
        "Binding_Position":   "{}-{}".format(bind_start, bind_end),
        "Mismatch_Count":     mm_count + gaps,
        "Similarity_Percent": similarity,
        "P_Value":            p_value,
        "All_Binding_Sites":  all_sites,
        "Fail_Reasons":       reasons,
        "Alignment_Structure": structure,
        "alignment":          alignment,
        # v3.0: Seed region MFE detayları (rule_seed_mfe aktifse dolu)
        "Seed_MFE_2_7":       seed_mfe_2_7,
        "Seed_MFE_8_13":      seed_mfe_8_13,
    }

if __name__ == "__main__":
    print("NanoSynapse Engine v2.1 -- diagnostics")
    mi = "UUAAUGCUAAUCGUGAUAGGGGU"
    # Perfect reverse complement of miRNA: ACCCCUAUCACGAUUAGCAUUAA
    tgt = "ACUGACACCCCUAUCACGAUUAGCAUUAACGUG"
    r = find_targets("hsa-miR-155-5p", mi, "test_gene", tgt, compute_pvalue=True, rule_set="animal")
    print("Animal Mode Status:   {}".format(r['Status']))
    print("Animal Confidence:    {}/100".format(r['Confidence_Score']))
    print("Animal MFE:           {} kcal/mol".format(r['MFE_kcal_mol']))
    print("Animal Penalty:       {}".format(r['Penalty_Score']))
    print("Animal Binding Pos:   {}".format(r['Binding_Position']))
    
    r_plant = find_targets("hsa-miR-155-5p", mi, "test_gene", tgt, compute_pvalue=True, rule_set="plant")
    print("Plant Mode Status:    {}".format(r_plant['Status']))
    print("Plant Confidence:     {}/100".format(r_plant['Confidence_Score']))
    print("Plant MFE:            {} kcal/mol".format(r_plant['MFE_kcal_mol']))
    print("Plant Penalty:        {}".format(r_plant['Penalty_Score']))
    print("Plant Binding Pos:    {}".format(r_plant['Binding_Position']))
    print("Engine ready [OK]")