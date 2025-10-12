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
        if (c === m_c) {
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
        const prematch = 'N'.repeat(m.index - start) + pattern_seq + 'N'.repeat(end - start - pattern_seq.length)
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
    let all_colliding_names = new Set()
    while (replacements.length > 0) {
        console.log(JSON.stringify(replacements))
        let new_replacements = []
        for (const replacement of replacements) {
            const seqs = replacement.seqs
            const replaced = replacement.replaced
            // pick exactly one index to change: try each possible next codon index (one per new branch)
            let id_start = 0
            if (replaced.length > 0)
                id_start = replaced[replaced.length - 1] + 1
            if (id_start >= ncodons)
                continue
            for (let i = id_start; i < ncodons; i++) {
                for (let s of seqs) {
                    let new_seqs = []
                    const codon = s.slice(i * 3, i * 3 + 3)
                    let aa = code.code[codon]
                    for (let new_codon of code.codons(aa)) {
                        if (new_codon == codon)
                            continue
                        let new_seq = s.slice(0, i * 3) + new_codon + s.slice(i * 3 + 3)
                        let matches = patterns_re.filter(p => new_seq.match(p.pattern))
                        if (matches.length == 0) {
                            if (n_mismatches == 1 || num_mismatches(new_seq, pattern) >= n_mismatches)
                                return {seq: new_seq, colliding_names: null};
                        }
                        let colliding_names = new Set(matches.map(c => c.name))
                        for (const n of colliding_names) all_colliding_names.add(n)
                        new_seqs.push(new_seq)
                    }
                    if (new_seqs.length > 0) {
                        new_replacements.push({seqs: new_seqs, replaced: [...replaced, i]})
                    }
                }
            }
            // let id_start = 0
            // if (replaced.length > 0)
            //     id_start = replaced[replaced.length - 1] + 1
            // if (id_start >= ncodons)
            //     continue
            // for (s of seqs) {
            //     for (let i = id_start; i < ncodons; i++) {
            //         let new_seqs = []
            //         const codon = s.slice(i * 3, i * 3 + 3)
            //         let aa = code.code[codon]
            //         for (new_codon of code.codons(aa)) {
            //             if (new_codon == codon)
            //                 continue
            //             let new_seq = s.slice(0, i * 3) + new_codon + s.slice(i * 3 + 3)                        
            //             let matches = patterns_re.filter(p => new_seq.match(p.pattern))
            //             if (matches.length == 0)
            //                 if (n_mismatches == 1 || num_mismatches(new_seq, pattern) >= n_mismatches)
            //                     return {seq: new_seq, colliding_names: null};
            //             let colliding_names = new Set(matches.map(c => c.name))
            //             all_colliding_names = all_colliding_names.union(colliding_names)
            //             new_seqs.push(new_seq)
            //         }
            //         if (new_seqs.length > 0)
            //             console.log(JSON.stringify([...replaced, i]))
            //             new_replacements.push({seqs: new_seqs, replaced: [...replaced, i]})
            //     }
            }
        // }
        replacements = new_replacements
    }
    return {seq: null, colliding_names: all_colliding_names}
}

function collision_patterns_re(collisions, code_left, code_right) {
    let res = new Set()
    function block_pattern(fit, rc=false) {
        return '.'.repeat(fit.start - code_left) + iupac_re(fit.pattern.seq) + '.'.repeat(code_right - (fit.end))
    }
    for (fit of collisions)
        res.add({pattern: block_pattern(fit), name: fit.pattern.name})
    return Array.from(res)
}

function try_remove_patterns(seq, patterns, code, n_mismatches=1) {
    const fits = try_insert_patterns(seq, patterns, code).slice().sort((f1, f2) => f1.start - f2.start)
    function find_colliding_fits(efit) {
        const coll = fits.filter(f => ((f.start >= efit.start) && (f.start < efit.end)) ||
                                      ((f.end >= efit.start)   && (f.end < efit.end)) ||
                                      ((efit.start >= f.start) && (efit.start < f.end)) ||
                                      ((efit.end >= f.start)   && (efit.end < f.end)))
        return Array.from(new Set(coll))
    }

    const exact_fits = fits.filter(f => f.exact)
    let history = []
    for (let efit of exact_fits) {
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
        
        let new_replacement = try_replace_codons(seq.slice(code_left, code_right), patterns_re, pattern, code, n_mismatches)
        let new_seq = new_replacement.seq
        if(new_seq !== null)
            seq = seq.slice(0, code_left) + new_seq + seq.slice(code_right)
        let record = {
            seq: seq,
            efit: efit,
            success: true,
            code_left: code_left,
            code_right: code_right,
            collisions: new_replacement.colliding_names
        }
        if (new_seq === null)
            record.success = false
        history.push(record)
        if (new_seq === null) 
            break
    }
    return history
}