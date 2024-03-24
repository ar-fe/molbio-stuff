function rev_comp(seq) {
    const comp = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N'
    }
    return seq.split('').reverse().map(c => comp[c]).join('')
}

const fwd_code = {
    'TTT':   'F',  'TAT':   'Y',  'TCT':   'S',  'TGT':   'C',  'ATT':   'I',  'AAT':   'N',  'ACT':   'T',  'AGT':   'S',
    'TTC':   'F',  'TAC':   'Y',  'TCC':   'S',  'TGC':   'C',  'ATC':   'I',  'AAC':   'N',  'ACC':   'T',  'AGC':   'S',
    'TTA':   'L',  'TAA':   '*',  'TCA':   'S',  'TGA':   '*',  'ATA':   'I',  'AAA':   'K',  'ACA':   'T',  'AGA':   'R',
    'TTG':   'L',  'TAG':   '*',  'TCG':   'S',  'TGG':   'W',  'ATG':   'M',  'AAG':   'K',  'ACG':   'T',  'AGG':   'R',
    'CTT':   'L',  'CAT':   'H',  'CCT':   'P',  'CGT':   'R',  'GTT':   'V',  'GAT':   'D',  'GCT':   'A',  'GGT':   'G',
    'CTC':   'L',  'CAC':   'H',  'CCC':   'P',  'CGC':   'R',  'GTC':   'V',  'GAC':   'D',  'GCC':   'A',  'GGC':   'G',
    'CTA':   'L',  'CAA':   'Q',  'CCA':   'P',  'CGA':   'R',  'GTA':   'V',  'GAA':   'E',  'GCA':   'A',  'GGA':   'G',
    'CTG':   'L',  'CAG':   'Q',  'CCG':   'P',  'CGG':   'R',  'GTG':   'V',  'GAG':   'E',  'GCG':   'A',  'GGG':   'G',
}

function translate(seq) {
    return (seq.match(/.{3}/g) || []).map(codon => fwd_code[codon]).join('')
}

const rev_code = {
    '*': ['TGA', 'TAA', 'TAG'], 'A': ['GCC', 'GCT', 'GCA', 'GCG'], 'M': ['ATG'], 'G': ['GGC', 'GGA', 'GGG', 'GGT'], 'S': ['AGC', 'TCC', 'TCT', 'TCA', 'AGT', 'TCG'],
    'C': ['TGC', 'TGT'], 'N': ['AAT', 'AAC'], 'H': ['CAC', 'CAT'], 'T': ['ACC', 'ACA', 'ACT', 'ACG'], 'D': ['GAC', 'GAT'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'I': ['ATC', 'ATT', 'ATA'], 'V': ['GTG', 'GTC', 'GTT', 'GTA'],'E': ['GAG', 'GAA'], 'Q': ['CAG', 'CAA'], 'K': ['AAG', 'AAA'], 'W': ['TGG'], 'F': ['TTC', 'TTT'],
    'R': ['AGA', 'AGG', 'CGG', 'CGC', 'CGA', 'CGT'], 'L': ['CTG', 'CTC', 'CTT', 'TTG', 'TTA', 'CTA'], 'Y': ['TAC', 'TAT']   		
}

function back_translate(seq) {
    return seq.split('').map(aa => rev_code[aa][0]).join('')
}

function prot_to_deg(seq) {
    const rev_code_deg = {
        '*': 'TRR',
        'A': 'GCN', 'M': 'ATG', 'G': 'GGN', 'S': 'WSN',
        'C': 'TGY', 'N': 'AAY', 'H': 'CAY', 'T': 'ACN',
        'D': 'GAY', 'P': 'CCN', 'I': 'ATH', 'V': 'GTN',
        'E': 'GAR', 'Q': 'CAR', 'K': 'AAR', 'W': 'TGG',
        'F': 'TTY', 'R': 'MGN', 'L': 'YTN', 'Y': 'TAY',
    }
    return seq.split('').map(a => rev_code_deg[a]).join('')
}

function pattern_to_prematch_re(pattern, add_rev_comp) {
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
    if (add_rev_comp) {
        const rc = rev_comp(pattern)
        if (rc != pattern) pattern_re += '|' + site_to_prematch(rc)
    }
    return '(?=(' + pattern_re + '))'
}

function try_fit_prematch(prematch, seq) {
    // prematch -- an original DNA seq with prematch place replaced by (possibly degenerate) RE site
    // seq -- the original seq. both arguments are minimal length extended to full codons
    const insert_patterns = {
        'A': '[A]', 'G': '[G]', 'C': '[C]', 'T': '[T]', 
        'R': '[AG]', 'Y': '[TC]', 'S': '[CG]', 'W': '[TA]',
        'K': '[TG]', 'M': '[AC]', 'B': '[TCG]', 'D': '[TAG]',
        'H': '[TAC]', 'V': '[ACG]', 'N': '[TACG]'
    }
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
        const m_c_pattern = insert_patterns[m_c[0]] + insert_patterns[m_c[1]] + insert_patterns[m_c[2]]
        if (c.match(m_c_pattern)) {
            result_codons.push(c)
            continue
        }
        for (const aa_c of rev_code[fwd_code[c]]) {
            if (aa_c.match(m_c_pattern)) {
                result_codons.push(aa_c)
                continue loop
            }
        }
        return null
    }
    return result_codons.join('')
}

function try_insert_pattern(seq, pattern, add_rev_comp) {
    // seq -- DNA seq, must be 3*n in length
    let fits = []
    const deg_seq = prot_to_deg(translate(seq))
    let prematches = deg_seq.matchAll(pattern_to_prematch_re(pattern, add_rev_comp))
    for (m of prematches) {
        const start = m.index - m.index % 3
        const end = m.index + pattern.length + 1 + (3 - (m.index - start)) 
        const subseq = seq.slice(start, end)
        const prematch = subseq.substring(0, m.index - start) + pattern + subseq.substring(m.index - start + pattern.length)
        let fit = try_fit_prematch(prematch, subseq)
        let site_rc = rev_comp(pattern)
        if (pattern != site_rc)
            if (!fit)
                fit = try_fit_prematch(subseq.substring(0, m.index - start) + site_rc + subseq.substring(m.index - start + pattern.length), subseq);
        if (fit) {
            fits.push({pos: m.index, seq: fit.substring(m.index - start, m.index - start + pattern.length)})
        }
    }
    return fits
}