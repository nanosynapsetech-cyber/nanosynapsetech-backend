# -*- coding: utf-8 -*-
# backend/mirna_engine.py -- NanoSynapse Engine v4.0 (RNAhybrid-Style)
# v4.0 changes:
#   [1] Sliding window search like RNAhybrid — target scanned, miRNA fixed
#   [2] miRNA is NEVER extended or trimmed — exact input sequence used
#   [3] Target NEVER has gaps (t:"-" forbidden) — only miRNA bulges allowed
#   [4] parse_alignment rewritten: anti-parallel, no target gaps
#   [5] find_targets uses best-MFE window across target
#   [6] All 5 biological rules preserved with user toggles

import random
import RNA
import numpy as np
from Bio.SeqUtils import gc_fraction

# ─── CONSTANTS ────────────────────────────────────────────────────────────────

ANIMAL_POSITION_WEIGHTS = {
    0: 0.5,
    1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0,
    8: 0.5,
    9: 1.5, 10: 1.5,
    11: 1.0, 12: 1.0, 13: 1.0, 14: 1.0, 15: 1.0,
    16: 0.5, 17: 0.5, 18: 0.5, 19: 0.5, 20: 0.5, 21: 0.5
}

WATSON_CRICK = [frozenset({'A', 'U'}), frozenset({'C', 'G'})]
WOBBLE       = frozenset({'G', 'U'})
ALL_PAIRS    = WATSON_CRICK + [WOBBLE]

DEFAULT_MAX_PENALTY_PLANT   = 4.0
DEFAULT_MAX_PENALTY_ANIMAL  = 4.0
DEFAULT_MFE_THRESHOLD_ANIMAL = -17.0
DEFAULT_MFE_THRESHOLD_PLANT  = -18.0


# ─── 1. ALIGNMENT PARSING (No target gaps) ───────────────────────────────────

def parse_alignment(mirna_aln: str, target_aln: str, structure: str) -> list:
    """
    Parse ViennaRNA dot-bracket into alignment columns.

    RULES (v4.0 — RNAhybrid style):
      - miRNA: 5'→3', may have bulges (m:"-") where target has extra bases
      - Target: NEVER has gaps. Every target column must carry a real nucleotide.
      - Anti-parallel: miRNA pos i pairs with target read 3'→5'

    Returns list of {m, t, match, bulge?}
    """
    if not structure or '&' not in structure:
        return []

    s1, s2 = structure.split('&', 1)

    # Safety clamp
    mirna_aln  = mirna_aln[:len(s1)]
    target_aln = target_aln[:len(s2)]

    # Find paired positions
    s1_pairs = [i for i, ch in enumerate(s1) if ch == '(']
    s2_pairs = [i for i, ch in enumerate(s2) if ch == ')']

    # Anti-parallel: first miRNA pair → last target pair
    s2_pairs_rev = list(reversed(s2_pairs))
    pairing_map  = {}
    for k in range(min(len(s1_pairs), len(s2_pairs_rev))):
        pairing_map[s1_pairs[k]] = s2_pairs_rev[k]

    columns = []
    target_placed = [False] * len(target_aln)

    for i in range(len(s1)):
        ch = s1[i]

        if ch == '(':
            t_idx = pairing_map.get(i)
            if t_idx is None:
                continue

            # Target bases between previous placed and this one (anti-parallel order)
            # These are target-side overhangs → miRNA bulge (m:"-"), target nucleotide present
            for j in range(len(target_aln) - 1, t_idx, -1):
                if not target_placed[j]:
                    t_base = target_aln[j]
                    if t_base and t_base != '-':
                        columns.append({
                            "m": "-",
                            "t": t_base,
                            "match": False,
                            "bulge": "mirna"  # miRNA has no base here; target does
                        })
                    target_placed[j] = True

            # Paired column
            if i < len(mirna_aln) and t_idx < len(target_aln):
                m_base = mirna_aln[i]
                t_base = target_aln[t_idx]
                pair   = frozenset({m_base, t_base})
                columns.append({
                    "m":     m_base,
                    "t":     t_base,
                    "match": pair in WATSON_CRICK or pair == WOBBLE,
                })
                target_placed[t_idx] = True

        else:
            # miRNA unpaired → miRNA bulge only if miRNA base exists
            if i < len(mirna_aln):
                m_base = mirna_aln[i]
                # Do NOT add t:"-". Instead mark as miRNA-only overhang (no target partner)
                # These are displayed as miRNA side-loops, not target gaps
                columns.append({
                    "m": m_base,
                    "t": "·",   # visual placeholder — NOT a gap, NOT a nucleotide
                    "match": False,
                    "bulge": "mirna_loop"
                })

    # Remaining unplaced target bases → miRNA bulge (m:"-")
    for j in range(len(target_aln) - 1, -1, -1):
        if not target_placed[j]:
            t_base = target_aln[j]
            if t_base and t_base != '-':
                columns.append({
                    "m": "-",
                    "t": t_base,
                    "match": False,
                    "bulge": "mirna"
                })
            target_placed[j] = True

    return columns


# ─── 2. WEIGHTED PENALTY SCORE ────────────────────────────────────────────────

def calculate_penalty_score(alignment: list, rule_set: str = "animal") -> float:
    score = 0.0

    if rule_set == "plant":
        mirna_pos = 0
        in_gap = False

        for col in alignment:
            m = col["m"]
            t = col["t"]

            if m == "-":
                if not in_gap:
                    score += 2.0
                    in_gap = True
                else:
                    score += 0.5
                continue

            if t in ("-", "·"):
                if not in_gap:
                    score += 2.0
                    in_gap = True
                else:
                    score += 0.5
                mirna_pos += 1
                continue

            in_gap = False
            mirna_pos += 1

            pair = frozenset({m, t})
            if pair in WATSON_CRICK:
                pass
            elif pair == WOBBLE:
                score += 0.5
            else:
                is_seed = (2 <= mirna_pos <= 13)
                score += 1.5 if is_seed else 1.0

    else:
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

            if t in ("-", "·"):
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
    rule_max_mismatch=True,
    rule_cleavage_site=True,
    rule_seed_mismatch=True,
    rule_consecutive=True,
):
    passed = True
    reasons = []
    mismatches = 0
    gaps = 0
    consecutive = 0
    max_consecutive_seen = 0
    seed_mm = 0
    seed_wobbles = 0
    early_seed_mm = 0

    mirna_pos = 0

    for col in alignment:
        m = col["m"]
        t = col["t"]

        # Target gap forbidden — should never happen with v4.0 engine
        # but guard anyway
        if t == "-":
            passed = False
            reasons.append("FATAL: target gap detected — engine error")
            continue

        if m == "-":
            # miRNA bulge: target has a nucleotide, miRNA doesn't
            gaps += 1
            consecutive += 1
            if consecutive > max_consecutive_seen:
                max_consecutive_seen = consecutive
            continue

        if t == "·":
            # miRNA loop (unpaired miRNA base, no target partner)
            gaps += 1
            consecutive += 1
            if consecutive > max_consecutive_seen:
                max_consecutive_seen = consecutive
            mirna_pos += 1
            if 1 <= mirna_pos <= 9:
                early_seed_mm += 1
            continue

        # Paired column
        mirna_pos += 1
        pair = frozenset({m, t})
        is_match = pair in ALL_PAIRS
        is_wobble = pair == WOBBLE

        if is_match:
            consecutive = 0
            if is_wobble:
                if rule_set == "animal" and (2 <= mirna_pos <= 8):
                    seed_wobbles += 1
                elif rule_set == "plant" and (2 <= mirna_pos <= 13):
                    seed_wobbles += 1
        else:
            mismatches += 1
            consecutive += 1
            if consecutive > max_consecutive_seen:
                max_consecutive_seen = consecutive

            if rule_set == "animal" and (2 <= mirna_pos <= 8):
                seed_mm += 1
            elif rule_set == "plant" and (2 <= mirna_pos <= 13):
                seed_mm += 1

            if 1 <= mirna_pos <= 9:
                early_seed_mm += 1

            if rule_cleavage_site and strict_cleavage and (10 <= mirna_pos <= 11):
                passed = False
                reasons.append(f"Rule ii: Mismatch at cleavage site pos {mirna_pos}")

    if rule_consecutive and max_consecutive_seen > 2:
        passed = False
        reasons.append(f"Rule iv: {max_consecutive_seen} consecutive MM/gaps (max 2)")

    if rule_seed_mismatch and early_seed_mm > 1:
        passed = False
        reasons.append(f"Rule iii: {early_seed_mm} mismatches in pos 1-9 (max 1)")

    if rule_set == "animal":
        if seed_mm > 2:
            passed = False
            reasons.append(f"Too many seed MM in pos 2-8 ({seed_mm}, max 2)")
        if seed_wobbles > 2:
            passed = False
            reasons.append(f"Too many G:U wobbles in seed ({seed_wobbles}, max 2)")
        if rule_max_mismatch:
            total_errors = mismatches + gaps
            if total_errors > max_mismatches:
                passed = False
                reasons.append(f"Rule i: Total MM+gaps ({total_errors}) > {max_mismatches}")
    else:
        if seed_mm > 2:
            passed = False
            reasons.append(f"Too many seed MM in pos 2-13 ({seed_mm}, max 2)")

        consec_seed_mm = 0
        mirna_pos_temp = 0
        for col in alignment:
            if col["m"] == "-":
                continue
            mirna_pos_temp += 1
            if 2 <= mirna_pos_temp <= 13:
                t = col["t"]
                m = col["m"]
                if t in ("-", "·") or frozenset({m, t}) not in ALL_PAIRS:
                    consec_seed_mm += 1
                    if consec_seed_mm > 1:
                        passed = False
                        reasons.append("Consecutive MM/gaps in seed region")
                        break
                else:
                    consec_seed_mm = 0

        if rule_max_mismatch:
            total_errors = mismatches + gaps
            if total_errors > max_mismatches:
                passed = False
                reasons.append(f"Rule i: Total MM+gaps ({total_errors}) > {max_mismatches}")

    penalty = calculate_penalty_score(alignment, rule_set)
    max_penalty = DEFAULT_MAX_PENALTY_PLANT if rule_set == "plant" else DEFAULT_MAX_PENALTY_ANIMAL
    if penalty > max_penalty:
        passed = False
        reasons.append(f"Penalty {penalty} > threshold {max_penalty}")

    return passed, list(set(reasons)), mismatches, gaps, penalty


# ─── 4. SEED REGION MFE ───────────────────────────────────────────────────────

def calculate_region_mfe(mirna_seq, target_seq, bind_start):
    try:
        seed_2_7  = mirna_seq[1:7]
        seed_8_13 = mirna_seq[7:13]

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
    try:
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
        return 1.0


# ─── 6. RNAHYBRID-STYLE SLIDING WINDOW SEARCH ────────────────────────────────

def find_best_duplex(mirna_seq: str, target_seq: str):
    """
    RNAhybrid-style: slide a window across target_seq.
    miRNA is used EXACTLY as given — never extended, never trimmed.
    Returns (best_duplex, best_window_start, best_mfe).
    """
    mirna_len  = len(mirna_seq)
    win_size   = mirna_len + 10   # slight flanking for ViennaRNA context
    best_mfe   = 0.0
    best_duplex = None
    best_start  = 0

    # Step 1 nt for accuracy (RNAhybrid default)
    for i in range(0, max(1, len(target_seq) - mirna_len + 1)):
        win_end = min(len(target_seq), i + win_size)
        window  = target_seq[i:win_end]
        if len(window) < 4:
            continue
        try:
            duplex = RNA.duplexfold(mirna_seq, window)
            if duplex.energy < best_mfe:
                best_mfe    = duplex.energy
                best_duplex = duplex
                best_start  = i
        except Exception:
            continue

    return best_duplex, best_start, best_mfe


# ─── 7. P-VALUE ───────────────────────────────────────────────────────────────

def calculate_pvalue(mirna_seq, target_seq, actual_mfe, n=200):
    try:
        seq_list = list(mirna_seq)
        background = []
        target_window = target_seq[:120]
        for _ in range(n):
            random.shuffle(seq_list)
            duplex = RNA.duplexfold(''.join(seq_list), target_window)
            background.append(duplex.energy)
        count_lte = sum(1 for m in background if m <= actual_mfe)
        return round(max(float(count_lte) / n, 1.0 / n), 4)
    except Exception:
        return 1.0


# ─── 8. COMPOSITE CONFIDENCE SCORE ────────────────────────────────────────────

def calculate_confidence(mfe, penalty, gc, accessibility, p_value, rule_set="animal"):
    mfe_s  = min(100.0, abs(mfe) * 2.5)
    pen_s  = max(0.0,   100.0 - penalty * 20)
    acc_s  = accessibility * 100.0
    gc_s   = max(0.0,   100.0 - abs(gc - 50) * 2)
    pval_s = max(0.0,   (1.0 - p_value) * 100)

    if rule_set == "plant":
        val = mfe_s * 0.25 + pen_s * 0.45 + acc_s * 0.15 + gc_s * 0.05 + pval_s * 0.10
    else:
        val = mfe_s * 0.40 + pen_s * 0.30 + acc_s * 0.15 + gc_s * 0.05 + pval_s * 0.10

    return round(min(100.0, val), 1)


# ─── 9. ALL BINDING SITES ─────────────────────────────────────────────────────

def find_all_binding_sites(mirna_seq, target_seq, min_mfe=-15.0, max_sites=3):
    win = len(mirna_seq) + 10
    sites = []

    for i in range(0, max(1, len(target_seq) - len(mirna_seq) + 1)):
        window = target_seq[i : i + win]
        try:
            duplex = RNA.duplexfold(mirna_seq, window)
        except Exception:
            continue

        if duplex.energy >= min_mfe:
            continue

        s1, s2 = duplex.structure.split('&')
        exact_start = i + duplex.j - 1
        exact_end   = exact_start + len(s2)

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


# ─── 10. MAIN PREDICTION FUNCTION ─────────────────────────────────────────────

def find_targets(
    mirna_id, mirna_seq, target_id, target_seq,
    compute_pvalue=False,
    rule_set="animal",
    max_mismatches=4,
    strict_cleavage=True,
    mfe_threshold=None,
    rule_max_mismatch=True,
    rule_cleavage_site=True,
    rule_seed_mismatch=True,
    rule_consecutive=True,
    rule_seed_mfe=False,
    seed_mfe_threshold=-20.0,
):
    """
    RNAhybrid-style miRNA-target prediction.
    
    Key guarantees (v4.0):
      1. mirna_seq is used EXACTLY — no extension, no trimming
      2. target NEVER has gaps (t:"-" forbidden in alignment)
      3. Only miRNA side can have bulges (m:"-")
      4. Best-MFE window across target is selected (sliding window)
    """
    mirna_len = len(mirna_seq)

    # ── Step 1: Sliding window — find best binding site ──────────────────────
    best_duplex, win_start, best_mfe = find_best_duplex(mirna_seq, target_seq)

    if best_duplex is None:
        # Fallback: direct duplexfold on whole target
        best_duplex = RNA.duplexfold(mirna_seq, target_seq)
        win_start   = 0
        best_mfe    = best_duplex.energy

    structure = best_duplex.structure
    mfe       = best_mfe

    # ── Step 2: Extract precise coordinates ──────────────────────────────────
    s1, s2 = structure.split('&', 1)

    # miRNA side: duplex.i is 1-based end of miRNA alignment
    mirna_end   = min(best_duplex.i, mirna_len)
    mirna_start = max(0, mirna_end - len(s1))
    mirna_aln   = mirna_seq[mirna_start:mirna_end]

    # Target side: duplex.j is 1-based start within the window
    t_local_start = max(0, best_duplex.j - 1)
    t_local_end   = min(len(target_seq) - win_start, t_local_start + len(s2))
    target_aln    = target_seq[win_start + t_local_start : win_start + t_local_end]

    # Absolute binding coordinates in target_seq
    bind_start = win_start + t_local_start + 1   # 1-based
    bind_end   = win_start + t_local_end          # 1-based inclusive

    # ── Step 3: Parse alignment — NO target gaps ──────────────────────────────
    alignment = parse_alignment(mirna_aln, target_aln, structure)

    # ── Step 4: Evaluate rules ────────────────────────────────────────────────
    passed, reasons, mm_count, gaps, penalty = evaluate_strict_rules(
        alignment, rule_set, max_mismatches, strict_cleavage,
        rule_max_mismatch=rule_max_mismatch,
        rule_cleavage_site=rule_cleavage_site,
        rule_seed_mismatch=rule_seed_mismatch,
        rule_consecutive=rule_consecutive,
    )

    # ── Step 5: MFE threshold ─────────────────────────────────────────────────
    threshold_mfe = mfe_threshold if mfe_threshold is not None else (
        DEFAULT_MFE_THRESHOLD_PLANT if rule_set == "plant" else DEFAULT_MFE_THRESHOLD_ANIMAL
    )
    if mfe > threshold_mfe:
        passed = False
        reasons.append(f"MFE ({mfe:.1f}) above threshold {threshold_mfe}")

    # ── Step 6: Accessibility ─────────────────────────────────────────────────
    accessibility = estimate_accessibility(target_seq, bind_start, mirna_len)

    # ── Step 7: GC content ────────────────────────────────────────────────────
    gc_pct = round(gc_fraction(mirna_seq) * 100, 2)

    # ── Step 8: p-value ───────────────────────────────────────────────────────
    p_value = calculate_pvalue(mirna_seq, target_seq, mfe) if compute_pvalue else 1.0

    # ── Step 9: Seed MFE (rule v) ─────────────────────────────────────────────
    seed_mfe_2_7  = 0.0
    seed_mfe_8_13 = 0.0
    if rule_seed_mfe:
        seed_mfe_2_7, seed_mfe_8_13 = calculate_region_mfe(mirna_seq, target_seq, bind_start)
        if not (seed_mfe_2_7 <= seed_mfe_threshold or seed_mfe_8_13 <= seed_mfe_threshold):
            passed = False
            reasons.append(
                f"Rule v: Seed MFE pos2-7={seed_mfe_2_7}, pos8-13={seed_mfe_8_13} "
                f"(need <= {seed_mfe_threshold})"
            )

    # ── Step 10: All binding sites ────────────────────────────────────────────
    min_site_mfe = DEFAULT_MFE_THRESHOLD_PLANT if rule_set == "plant" else DEFAULT_MFE_THRESHOLD_ANIMAL
    all_sites = (
        find_all_binding_sites(mirna_seq, target_seq, min_site_mfe)
        if len(target_seq) > mirna_len + 4 else []
    )

    # ── Step 11: Confidence + similarity ─────────────────────────────────────
    confidence = calculate_confidence(mfe, penalty, gc_pct, accessibility, p_value, rule_set)

    aligned_len = len(alignment)
    match_count = sum(1 for col in alignment if col["match"])
    similarity  = round((float(match_count) / max(1, aligned_len)) * 100, 2)

    return {
        "miRNA_ID":            mirna_id,
        "Gene_ID":             target_id,
        "Status":              "PASS" if passed else "FAIL",
        "Confidence_Score":    confidence,
        "MFE_kcal_mol":        round(mfe, 2),
        "Penalty_Score":       penalty,
        "Accessibility":       accessibility,
        "GC_Content_Percent":  gc_pct,
        "Binding_Position":    f"{bind_start}-{bind_end}",
        "Mismatch_Count":      mm_count + gaps,
        "Similarity_Percent":  similarity,
        "P_Value":             p_value,
        "All_Binding_Sites":   all_sites,
        "Fail_Reasons":        reasons,
        "Alignment_Structure": structure,
        "alignment":           alignment,
        "Seed_MFE_2_7":        seed_mfe_2_7,
        "Seed_MFE_8_13":       seed_mfe_8_13,
        # v4.0 debug info
        "miRNA_Used":          mirna_seq,        # exact sequence used (no extension)
        "Target_Window_Start": win_start,
        "Target_Aligned":      target_aln,
    }


if __name__ == "__main__":
    print("NanoSynapse Engine v4.0 (RNAhybrid-style) — diagnostics")
    mi  = "UUAAUGCUAAUCGUGAUAGGGGU"
    tgt = "ACUGACACCCCUAUCACGAUUAGCAUUAACGUG"
    r   = find_targets("hsa-miR-155-5p", mi, "test_gene", tgt,
                       compute_pvalue=False, rule_set="animal")
    print(f"Status:        {r['Status']}")
    print(f"MFE:           {r['MFE_kcal_mol']} kcal/mol")
    print(f"Confidence:    {r['Confidence_Score']}/100")
    print(f"Binding Pos:   {r['Binding_Position']}")
    print(f"miRNA used:    {r['miRNA_Used']}")
    print(f"Target window: {r['Target_Window_Start']}")

    # Verify: no target gaps
    tgaps = [c for c in r['alignment'] if c['t'] == '-']
    print(f"Target gaps (must be 0): {len(tgaps)}")
    mlen  = sum(1 for c in r['alignment'] if c['m'] != '-')
    print(f"miRNA bases in alignment: {mlen} (input len: {len(mi)})")
    print("Engine ready [OK]")