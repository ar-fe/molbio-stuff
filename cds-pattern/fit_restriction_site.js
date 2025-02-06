class Pattern {
    constructor(name, pattern_seq) {
        this.name = name
        this.seq = pattern_seq
    }
}

class Fit {
    constructor(pattern, pos, new_seq, exact) {
        this.pattern = pattern
        this.pos = pos
        this.start = pos
        this.new_seq = new_seq
        this.end = pos + new_seq.length
        this.exact = exact
    }
}

class GeneticCode {
    constructor(transl_table, freqs) {
        this.code = transl_table
        this.rev_code = Object.fromEntries(Array.from(new Set(Object.values(this.code))).map(aa => [aa, Object.keys(this.code).filter(codon => this.code[codon] == aa)]))
        if (freqs !== undefined)
            Object.values(this.rev_code).forEach(codons => codons.sort((c1, c2) => freqs[c2] - freqs[c1]))
        function codons_to_deg(codons) {
            const deg_map = {
                'A': 'A', 'G': 'G', 'C': 'C', 'T': 'T', 'AG': 'R', 'CT': 'Y', 'CG': 'S', 'AT': 'W', 'GT': 'K', 'AC': 'M',
                'CGT': 'B', 'AGT': 'D', 'ACT': 'H', 'ACG': 'V', 'ACGT': 'N',
            }
            let sets = [new Set(), new Set(), new Set()]
            for (let codon of codons)
                for (let i = 0; i < 3; i++)
                    sets[i].add(codon[i])
            return sets.map(set => deg_map[Array.from(set).sort().join('')]).join('')
        }
        this.deg_code = Object.fromEntries(Object.keys(this.rev_code).map(aa => [aa, codons_to_deg(this.rev_code[aa])]))
    }
    translate(cds) {
        return (cds.match(/.{3}/g) || []).map(codon => this.code[codon]).join('')
    }
    codons(aa) {
        return this.rev_code[aa]
    }
    back_translate(seq) {
        return seq.split('').map(aa => this.rev_code[aa][0]).join('')
    }
    prot_to_deg(prot_seq) {
        return prot_seq.split('').map(a => this.deg_code[a]).join('')
    }
    codon_variants(codon) {
        return this.rev_code[this.code[codon]]
    }
}

function rev_comp(seq) {
    const comp = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'
    }
    return seq.split('').reverse().map(c => comp[c]).join('')
}

// const fwd_code = {
//     'TTT':   'F',  'TAT':   'Y',  'TCT':   'S',  'TGT':   'C',  'ATT':   'I',  'AAT':   'N',  'ACT':   'T',  'AGT':   'S',
//     'TTC':   'F',  'TAC':   'Y',  'TCC':   'S',  'TGC':   'C',  'ATC':   'I',  'AAC':   'N',  'ACC':   'T',  'AGC':   'S',
//     'TTA':   'L',  'TAA':   '*',  'TCA':   'S',  'TGA':   '*',  'ATA':   'I',  'AAA':   'K',  'ACA':   'T',  'AGA':   'R',
//     'TTG':   'L',  'TAG':   '*',  'TCG':   'S',  'TGG':   'W',  'ATG':   'M',  'AAG':   'K',  'ACG':   'T',  'AGG':   'R',
//     'CTT':   'L',  'CAT':   'H',  'CCT':   'P',  'CGT':   'R',  'GTT':   'V',  'GAT':   'D',  'GCT':   'A',  'GGT':   'G',
//     'CTC':   'L',  'CAC':   'H',  'CCC':   'P',  'CGC':   'R',  'GTC':   'V',  'GAC':   'D',  'GCC':   'A',  'GGC':   'G',
//     'CTA':   'L',  'CAA':   'Q',  'CCA':   'P',  'CGA':   'R',  'GTA':   'V',  'GAA':   'E',  'GCA':   'A',  'GGA':   'G',
//     'CTG':   'L',  'CAG':   'Q',  'CCG':   'P',  'CGG':   'R',  'GTG':   'V',  'GAG':   'E',  'GCG':   'A',  'GGG':   'G',
// }

// function translate(seq) {
//     return (seq.match(/.{3}/g) || []).map(codon => fwd_code[codon]).join('')
// }

const human_freqs = {
    'TTT': 1385301,
    'TTC': 1413268,
    'TTA': 703680,
    'TTG': 1086777,
    'CTT': 1138433,
    'CTC': 1439345,
    'CTA': 601662,
    'CTG': 2918400,
    'ATT': 1331901,
    'ATC': 1508988,
    'ATA': 652939,
    'ATG': 1739992,
    'GTT': 949137,
    'GTC': 1086717,
    'GTA': 618960,
    'GTG': 2090923,
    'TAT': 978774,
    'TAC': 1090514,
    'TAA': 35218,
    'TAG': 28499,
    'CAT': 956479,
    'CAC': 1184041,
    'CAA': 1136523,
    'CAG': 2872161,
    'AAT': 1489775,
    'AAC': 1478832,
    'AAA': 2221062,
    'AAG': 2567940,
    'GAT': 1942185,
    'GAC': 1961667,
    'GAA': 2719693,
    'GAG': 3206546,
    'TCT': 1368632,
    'TCC': 1399962,
    'TCA': 1142684,
    'TCG': 325925,
    'CCT': 1560898,
    'CCC': 1544626,
    'CCA': 1529004,
    'CCG': 503096,
    'ACT': 1152700,
    'ACC': 1442511,
    'ACA': 1335468,
    'ACG': 452037,
    'GCT': 1534685,
    'GCC': 2088762,
    'GCA': 1377145,
    'GCG': 477758,
    'TGT': 841042,
    'TGC': 873765,
    'TGA': 63801,
    'TGG': 937286,
    'CGT': 367659,
    'CGC': 704401,
    'CGA': 518818,
    'CGG': 871786,
    'AGT': 1135376,
    'AGC': 1591829,
    'AGA': 1073213,
    'AGG': 980476,
    'GGT': 875715,
    'GGC': 1599325,
    'GGA': 1384137,
    'GGG': 1240793,
}


// function code_to_rev(code, freqs) {
//     let rev_code = Object.fromEntries(Array.from(new Set(Object.values(code))).map(aa => [aa, Object.keys(code).filter(codon => code[codon] == aa)]))
//     if (freqs !== undefined) {
//         Object.values(rev_code).forEach(codons => codons.sort((c1, c2) => freqs[c2] - freqs[c1]))
//     }
//     return rev_code
// }

// const rev_code = {
//     '*': ['TGA', 'TAA', 'TAG'], 'A': ['GCC', 'GCT', 'GCA', 'GCG'], 'M': ['ATG'], 'G': ['GGC', 'GGA', 'GGG', 'GGT'], 'S': ['AGC', 'TCC', 'TCT', 'TCA', 'AGT', 'TCG'],
//     'C': ['TGC', 'TGT'], 'N': ['AAT', 'AAC'], 'H': ['CAC', 'CAT'], 'T': ['ACC', 'ACA', 'ACT', 'ACG'], 'D': ['GAC', 'GAT'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
//     'I': ['ATC', 'ATT', 'ATA'], 'V': ['GTG', 'GTC', 'GTT', 'GTA'], 'E': ['GAG', 'GAA'], 'Q': ['CAG', 'CAA'], 'K': ['AAG', 'AAA'], 'W': ['TGG'], 'F': ['TTC', 'TTT'],
//     'R': ['AGA', 'AGG', 'CGG', 'CGC', 'CGA', 'CGT'], 'L': ['CTG', 'CTC', 'CTT', 'TTG', 'TTA', 'CTA'], 'Y': ['TAC', 'TAT']   		
// }

// function back_translate(seq) {
//     return seq.split('').map(aa => rev_code[aa][0]).join('')
// }

// function prot_to_deg(seq) {
//     const rev_code_deg = {
//         '*': 'TRR',
//         'A': 'GCN', 'M': 'ATG', 'G': 'GGN', 'S': 'WSN',
//         'C': 'TGY', 'N': 'AAY', 'H': 'CAY', 'T': 'ACN',
//         'D': 'GAY', 'P': 'CCN', 'I': 'ATH', 'V': 'GTN',
//         'E': 'GAR', 'Q': 'CAR', 'K': 'AAR', 'W': 'TGG',
//         'F': 'TTY', 'R': 'MGN', 'L': 'YTN', 'Y': 'TAY',
//     }
//     return seq.split('').map(a => rev_code_deg[a]).join('')
// }

function pattern_to_prematch_re(pattern) {
    const prematch_patterns = {
        'A': '[ARWMDHVN]', 'G': '[GRSKBDVN]', 'C': '[CYSMBHVN]', 'T': '[TYWKBDHN]',
        'R': '[AGRSWKMBDHVN]', 'Y': '[CTYSWKMBDHVN]', 'S': '[GCRYSKMBDHVN]', 'W': '[ATRYWKMBDHVN]', 
        'K': '[GTRYSWKBDHVN]', 'M': '[ACRYSWMBDHVN]', 'B': '[GCTRYSWKMBDHVN]', 'D': '[AGTRYSWKMBDHVN]', 
        'H': '[ACTRYSWKMBDHVN]', 'V': '[AGCRYSWKMBDHVN]', 'N': '[AGCTRYSWKMBDHVN]'
    }    
    function site_to_prematch(seq) {
        return seq.split('').map(nt => prematch_patterns[nt]).join('')
    }
    let pattern_re = site_to_prematch(pattern)
    return '(?=(' + pattern_re + '))'
}

function iupac_re(seq) {
    const iupac_re_map = {
        'A': 'A', 'G': 'G', 'C': 'C', 'T': 'T',
        'R': '[AG]', 'Y': '[TC]', 'S': '[CG]', 'W': '[TA]',
        'K': '[TG]', 'M': '[AC]', 'B': '[TCG]', 'D': '[TAG]',
        'H': '[TAC]', 'V': '[ACG]', 'N': '.'
    }
    return seq.split('').map(n => iupac_re_map[n]).join('')
}

function try_fit_prematch(prematch, seq, code) {
    // prematch -- an original DNA seq with prematch place replaced by (possibly degenerate) RE site
    // seq -- the original seq. both arguments are minimal length extended to full codons
    const m_codons = prematch.match(/.{3}/g)
    const codons = seq.match(/.{3}/g)
    let result_codons = []
    loop: for (let i = 0; i < codons.length; i++) {
        const m_c = m_codons[i]
        const c = codons[i]
        if (c == m_c) {
            result_codons.push(c)
            continue
        }
        const m_c_pattern = iupac_re(m_c)
        if (c.match(m_c_pattern)) {
            result_codons.push(c)
            continue
        }
        for (const aa_c of code.codon_variants(c)) {
            if (aa_c.match(m_c_pattern)) {
                result_codons.push(aa_c)
                continue loop
            }
        }
        return null
    }
    return result_codons.join('')
}

function try_insert_pattern(seq, deg_seq, pattern, code) {
    // seq -- DNA seq, must be 3*n in length
    let pattern_seq = pattern.seq
    let fits = []
    let prematches = deg_seq.matchAll(pattern_to_prematch_re(pattern_seq))
    for (m of prematches) {
        const start = m.index - m.index % 3
        const end = m.index + pattern_seq.length + 1 + (3 - (m.index - start)) 
        const subseq = seq.slice(start, end)
        const prematch = subseq.substring(0, m.index - start) + pattern_seq + subseq.substring(m.index - start + pattern_seq.length)
        let fit = try_fit_prematch(prematch, subseq, code)
        if (fit)
            fits.push(new Fit(pattern, m.index, fit.substring(m.index - start, m.index - start + pattern_seq.length), !!subseq.match(iupac_re(pattern_seq))))
    }
    return fits
}

function try_insert_patterns(seq, patterns, code) {
    let found_pattern_seqs = []
    let found_fits = []
    const deg_seq = code.prot_to_deg(code.translate(seq))
    for (const pattern of patterns) {
        if (found_pattern_seqs.includes(pattern.seq))
            continue
        fits = try_insert_pattern(seq, deg_seq, pattern, code)    
        if (fits.length > 0) {
            found_fits = found_fits.concat(fits)
            found_pattern_seqs.push(pattern.seq)
        }
    }
    return found_fits
}



function bin_search_fit_end(fits, value) {
    function bin_search(start, stop) {
        if (stop - start == 1) {
            if (fits[start].end <= value)
                return stop
            return start
        }
        let pivot = start + Math.floor((stop - start) / 2)
        if (fits[pivot].end <= value)
            return bin_search(pivot, stop)
        else
            return bin_search(0, pivot)
    }
    return bin_search(0, fits.length)
}

function try_replace_codons(seq, patterns_re, pattern, code, n_mismatches=1) {
    function num_mismatches(seq, pattern) {        
        const nuc_mismatches = {
            'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG', 'R': 'CT', 'Y': 'AG', 'S': 'AT', 'W': 'CG',
            'K': 'AC', 'M': 'GT', 'B': 'A', 'D': 'C', 'H': 'G', 'V': 'T', 'N': ''
        }
        let n_mm = 0
        for (let i = 0; i < seq.length; i++)
            if (nuc_mismatches[pattern[i]].includes(seq[i]))
                n_mm++;
        return n_mm
    }
    let ncodons = seq.length / 3
    let replacements = [{seqs: [seq], replaced: []}]
    while (replacements.length > 0) {
        console.log(replacements)
        let new_replacements = []
        for (const replacement of replacements) {
            const seqs = replacement.seqs
            const replaced = replacement.replaced
            let id_start = 0
            if (replaced.length > 0)
                id_start = replaced[replaced.length - 1] + 1
            if (id_start >= ncodons)
                continue
            for (seq of seqs) {
                for (let i = id_start; i < ncodons; i++) {
                    let new_seqs = []
                    const codon = seq.slice(i * 3, i * 3 + 3)
                    let aa = code.code[codon]
                    for (new_codon of code.codons(aa)) {
                        if (new_codon == codon)
                            continue
                        let new_seq = seq.slice(0, i * 3) + new_codon + seq.slice(i * 3 + 3)
                        if (!new_seq.match(patterns_re))
                            if (n_mismatches == 1 || num_mismatches(new_seq, pattern) >= n_mismatches)
                                return new_seq;
                        new_seqs.push(new_seq)
                    }
                    if (new_seqs.length > 0)
                        new_replacements.push({seqs: new_seqs, replaced: [...replaced, i]})
                }
            }
        }
        replacements = new_replacements
    }
    return null
}

function collision_patterns_re(collisions, code_left, code_right) {
    let res = new Set()
    function block_pattern(fit, rc=false) {
        return '.'.repeat(fit.start - code_left) + iupac_re(fit.pattern.seq) + '.'.repeat(code_right - (fit.end))
    }
    for (fit of collisions)
        res.add(block_pattern(fit))
    return Array.from(res).join('|')
}

function try_remove_patterns(seq, patterns, code, n_mismatches=1) {
    const fits = try_insert_patterns(seq, patterns, code)
    const fits_bystart = fits.slice().sort((f1, f2) => f1.start - f2.start)
    const fits_byend = fits.slice().sort((f1, f2) => f1.end - f2.end)
    function find_colliding_fits(efit) {
        const ends_left  = bin_search_fit_end(fits_byend, efit.start),
              ends_right = bin_search_fit_end(fits_byend, efit.end)
        return Array.from(new Set(fits_byend.slice(ends_left, ends_right)))
    }

    const exact_fits = fits_bystart.filter(f => f.exact)
    let history = []
    for (efit of exact_fits) {
        // if the sequence is already destroyed
        if (!seq.substring(efit.start, efit.end).match(iupac_re(efit.pattern.seq)))
            continue
        const collisions = find_colliding_fits(efit)
        const block_left = Math.min(efit.start, ...collisions.map(f => f.start)),
              block_right = Math.max(efit.end, ...collisions.map(f => f.end))
        // full codons containing the collision block
        const code_left = block_left - (block_left % 3),
              code_right = block_right + ((block_right % 3 == 0) ? 0 : (3 - block_right % 3))
        const patterns_re = collision_patterns_re(collisions, code_left, code_right)
        const pattern = 'N'.repeat(efit.start - code_left) + efit.pattern.seq + 'N'.repeat(code_right - efit.end)
        
        let new_seq = try_replace_codons(seq.slice(code_left, code_right), patterns_re, pattern, code, n_mismatches)

        if(new_seq !== null)
            seq = seq.slice(0, code_left) + new_seq + seq.slice(code_right)
        let record = {seq: seq, efit: efit, success: true, code_left: code_left, code_right: code_right} //, __dbg_p: patterns_re}
        if (new_seq === null)
            record.success = false
        history.push(record)
        if (new_seq === null) 
            break
    }
    return history
}