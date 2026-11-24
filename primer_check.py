# Primer QC: nearest-neighbour melting temperature, GC clamp, hairpins and self-dimers.
import math

# SantaLucia (1998) unified nearest-neighbour parameters: dH (kcal/mol), dS (cal/K/mol).
NN = {
    "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3), "CA": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "CT": (-7.8, -21.0), "GA": (-8.2, -22.2), "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4), "GG": (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3, 4.1)
R = 1.987
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def revcomp(seq):
    return "".join(COMP[b] for b in reversed(seq))


def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)


def tm_wallace(seq):
    # The "2 + 4" rule of thumb: fine for 14-mers, increasingly wrong after that.
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return 2 * at + 4 * gc


def tm_nearest_neighbour(seq, na_mM=50.0, primer_nM=250.0):
    seq = seq.upper()
    dh, ds = 0.0, 0.0
    for end in (seq[0], seq[-1]):
        h, s = INIT_GC if end in "GC" else INIT_AT
        dh, ds = dh + h, ds + s
    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        h, s = NN[pair] if pair in NN else NN[revcomp(pair)]
        dh, ds = dh + h, ds + s
    ds += 0.368 * (len(seq) - 1) * math.log(na_mM / 1000.0)
    ct = primer_nM * 1e-9
    return 1000.0 * dh / (ds + R * math.log(ct / 4.0)) - 273.15


def longest_homopolymer(seq):
    best, run = 1, 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def longest_hairpin_stem(seq, min_loop=3):
    # Longest run of i..i+k pairing with j..j-k, leaving a loop of at least min_loop.
    best = 0
    n = len(seq)
    for i in range(n):
        for j in range(n - 1, i + min_loop, -1):
            k = 0
            while (i + k < j - k - min_loop and COMP[seq[i + k]] == seq[j - k]):
                k += 1
            best = max(best, k)
    return best


def three_prime_dimer(a, b, window=5):
    # Longest contiguous complementary run that involves the last `window` bases of a.
    rb = b[::-1]
    best = 0
    for shift in range(-len(a) + 1, len(rb)):
        run = 0
        for i in range(len(a)):
            j = i + shift
            if 0 <= j < len(rb) and COMP[a[i]] == rb[j]:
                run += 1
                if i >= len(a) - window:
                    best = max(best, run)
            else:
                run = 0
    return best


def check_primer(seq, na_mM=50.0, primer_nM=250.0):
    seq = seq.upper()
    clamp = sum(b in "GC" for b in seq[-5:])
    report = {
        "primer": seq,
        "length": len(seq),
        "gc": round(gc_content(seq), 3),
        "tm_wallace": tm_wallace(seq),
        "tm_nn": round(tm_nearest_neighbour(seq, na_mM, primer_nM), 1),
        "gc_clamp": clamp,
        "homopolymer": longest_homopolymer(seq),
        "hairpin_stem": longest_hairpin_stem(seq),
        "self_dimer_3p": three_prime_dimer(seq, seq),
    }
    flags = []
    if not 18 <= report["length"] <= 25:
        flags.append("length")
    if not 0.40 <= report["gc"] <= 0.60:
        flags.append("gc")
    if not 52 <= report["tm_nn"] <= 65:
        flags.append("tm")
    if not 1 <= clamp <= 3:
        flags.append("clamp")
    if report["homopolymer"] >= 5:
        flags.append("homopolymer")
    if report["hairpin_stem"] >= 4:
        flags.append("hairpin")
    if report["self_dimer_3p"] >= 4:
        flags.append("3p_dimer")
    report["flags"] = ";".join(flags)
    report["pass"] = not flags
    return report


def check_pair(fwd, rev, max_tm_diff=3.0):
    a, b = check_primer(fwd), check_primer(rev)
    return {
        "tm_diff": round(abs(a["tm_nn"] - b["tm_nn"]), 1),
        "cross_dimer_3p": max(three_prime_dimer(fwd, rev), three_prime_dimer(rev, fwd)),
        "pair_ok": a["pass"] and b["pass"] and abs(a["tm_nn"] - b["tm_nn"]) <= max_tm_diff,
    }
