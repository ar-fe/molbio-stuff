
function make_codes_table() {
    let select = document.querySelector('#gen-code')
    let optgroup = select.querySelector('optgroup')

    for (taxon of codon_freqs) {
        let {taxid, name, freqs} = taxon
        optgroup.innerHTML += `<option data-transl-table="1" data-taxid="${taxid}">1. Standard (${name})</option>`
    }
    for (id in transl_tables) {
        if (id == 1) continue;
        select.innerHTML += `<option data-transl-table="${id}">${id}. ${transl_table_names[id]}</option>`
    }
}

make_codes_table()

function table_header(seq, translated) {
    const clustalX_colorscheme = {
        '*': '#ffffff',
        'A': '#80a0f0', 'D': '#c048c0', 'W': '#80a0f0', 'C': '#f08080',
        'I': '#80a0f0', 'N': '#15c015', 'V': '#80a0f0', 'G': '#f09048',
        'L': '#80a0f0', 'Q': '#15c015', 'K': '#f01505', 'P': '#c0c000',
        'M': '#80a0f0', 'S': '#15c015', 'R': '#f01505', 'H': '#15a4a4',
        'F': '#80a0f0', 'T': '#15c015', 'E': '#c048c0', 'Y': '#15a4a4',
    }
    let numbers = ''
    let start_pos = 1
    for (let i = start_pos; i <= translated.length + start_pos - 1; i += 5) {
        let s = String(i)
        if (i % 10 == 5 || numbers == '') {
            if (numbers == '') {
                numbers = s
                i -= i % 5
            } else {
                const to5 = '                . '
                numbers += to5.substring(Math.max(to5.length - (i - start_pos + 1) * 3 + s.length % 15, s.length % 15))
            }            
        } else {
            numbers += '               '.substring(0, Math.min((i - start_pos) * 3 - numbers.length % 15, 27 - numbers.length % 30))
            if(numbers[numbers.length - 1] == ' ') numbers += s
        }        
    }
    numbers = `<span style="color: gray; user-select: none;">${numbers}</span>`
    let prot_seq = translated.split('')
        .map(c => `<span style="background-color: ${clustalX_colorscheme[c]};"> ` + c + ' </span>').join('')
    let tr = document.createElement('tr')
    tr.innerHTML = '<th id="table-topleft">Click<br/>to →scroll:</th>' + 
        '<th>' + numbers + '<br/>' + prot_seq + '<br/>' + seq + '</th>'
    return tr
}

function get_input_filtered(id, filter) {
    let input = document.getElementById('in-seq').value.toUpperCase()
    return input.split('').filter(c => filter.includes(c)).join('')
}

function get_input_cds(code) {
    let input_type = document.querySelector('input[name="input-type"]:checked').id
    let cds = ''
    if (input_type == 'in-prot-radio') {
        const prot = get_input_filtered('in-seq', 'ACDEFGHIKLMNPQRSTVWY')
        cds = code.back_translate(prot)
    } else {
        cds = get_input_filtered('in-seq', 'AGTC')
        cds = cds.substring(0, cds.length - cds.length % 3)
        if (document.getElementById('codon-optimize').checked)
            cds = code.back_translate(code.translate(cds))
    }
    return cds
}

function get_input_patterns() {
    function filter_seq(seq) {
        return seq.toUpperCase().split('').filter(c => 'AGCTRYSWKMBDHVN'.includes(c)).join('')
    }
    function add_rc_patterns(patterns) {
        return patterns.concat(patterns.filter(p => p.seq != rev_comp(p.seq)).map(p => new Pattern(p.name, rev_comp(p.seq)))).sort((p1, p2) => p1.name.localeCompare(p2.name))
    }
    let patterns_custom = [], patterns_re = []
    if (document.querySelector('#custom-seq-checkbox').checked) {
        patterns_custom = document.querySelector('#custom-seq')
            .value.split('\n').filter(l => l != '').map(s => s.split(':'))
            .map(np => np.length == 2 ? new Pattern(np[0], filter_seq(np[1])) : new Pattern(filter_seq(np[0]), filter_seq(np[0])))
        if (document.querySelector('#seqs-bottom-strand').checked)
            patterns_custom = add_rc_patterns(patterns_custom)
    }
    if (document.querySelector('#re-checkbox').checked) {
        patterns_re = Array.from(document.querySelectorAll('#restriction-enzymes input:checked')).map(el => new Pattern(el.id, restriction_enzymes[el.id]))
        patterns_re = add_rc_patterns(patterns_re)
    }
    return [...patterns_re, ...patterns_custom]
}

function highlight_site_insert(seq, fit) {
    return seq.substring(0, fit.pos) +
        '<span class="highlight highlight-insert">' +
        fit.new_seq +
        '</span>' +
        seq.substring(fit.end)
}

function highlight_site_remove(seq, record) {
    return record.seq.substring(0, record.code_left) +
        '<span class="highlight highlight-remove-outer">' +
        record.seq.substring(record.code_left, record.efit.start) +
        '<span class="highlight highlight-remove">' +
        record.seq.substring(record.efit.start, record.efit.end) +
        '</span>' +
        record.seq.substring(record.efit.end, record.code_right) +
        '</span>' +
        record.seq.substring(record.code_right)
}

function result_tr(fit, highlighted) {
    let tr = document.createElement('tr')
    tr.innerHTML = 
        `<td><span class="action-link" onclick="scroll_to_match(this)" title="${fit.pattern.seq}">${fit.pattern.name}</span></td>` +
        `<td>${highlighted}</td>`
    return tr
}
function add_status_string(str) {
    let status = document.createElement('div')
    status.classList.add('status-string')
    status.innerHTML = `<td colspan="2">${str}</td>`
    document.getElementById('results').appendChild(status)
}

function run_insert(seq, patterns, code) {
    let fits = try_insert_patterns(seq, patterns, code)
    fits.sort((f1, f2) => ((f1.pos - f2.pos) == 0) ? (f1.pattern.name.localeCompare(f2.pattern.name)) : (f1.pos - f2.pos))
    fits.map(fit => result_tr(fit, highlight_site_insert(seq, fit))).forEach(tr => document.querySelector('table').appendChild(tr))
    if (fits.length == 0) {
        add_status_string('<span style="color: red;">Nothing found :(</span>')
    }
}

function run_remove(seq, patterns, code) {
    let history = try_remove_patterns(seq, patterns, code, +document.querySelector('#n-replacements').value, code)
    history.map(record => result_tr(record.efit, highlight_site_remove(seq, record))).forEach(tr => document.querySelector('table').appendChild(tr))
    let last = history[history.length - 1]
    if (history.length > 0 && last.success == true)
        add_status_string(`<span style="color: green;">Succesfully removed ${history.length} sequence(s)</span>`)
    else if(history.length == 0)
        add_status_string(`<span style="color: green;">Nothing to remove</span>`)
    else {
        let len = Array.from(last.collisions).length
        if (len === 1 && last.collisions.has(last.efit.pattern.name))
            add_status_string(`<span style="color: red;">Unable to remove "${last.efit.pattern.name}"</span>`)
        else {
            last.collisions.delete(last.efit.pattern.name)
            let extra_collisions = Array.from(last.collisions)
            add_status_string(`<span style="color: red;">Unable to remove "${last.efit.pattern.name}"` + 
                ` (collision(s) with ${extra_collisions.join(', ')})</span>`)
        }
    }
}
function get_code() {
    const selected_code = document.querySelector('#gen-code').options[document.querySelector('#gen-code').selectedIndex]
    return new GeneticCode(transl_tables[selected_code.dataset.translTable], codon_freqs.filter(tax => tax.taxid == selected_code.dataset.taxid)[0].freqs)
}

function calculate() {
    const code = get_code()
    let cds = get_input_cds(code)
    let patterns = get_input_patterns()
    
    let results_table = document.getElementById('results')
    results_table.innerHTML = ''
    if (cds == '')
        return
    results_table.appendChild(table_header(cds, code.translate(cds)))
    if (patterns.length == 0)
        return
    if (document.querySelector('input[name=mode]:checked').id == 'mode-insert')
        run_insert(cds, patterns, code)
    else
        run_remove(cds, patterns, code)
    results_table.scrollLeft = 0
    results_table.scrollTop = 0
}

function rev_comp_input() {
    document.getElementById('in-seq').value = rev_comp(get_input_cds(get_code()))
}

let table = document.getElementById('restriction-enzymes')
for (RE of Object.keys(restriction_enzymes)) {
    let tr = document.createElement('tr')
    tr.innerHTML = `<td><input type="checkbox" name="re-list" id="${RE}" ${subsets.sba.includes(RE) ? "checked" : ""}><label for="${RE}">${RE}</label></td><td style="font-family: monospace;">${restriction_enzymes[RE]}</td>`
    table.appendChild(tr)
}

function scroll_to_match(el) {
    let highlight = el.parentNode.parentNode.querySelector('.highlight')
    document.getElementById('results').scrollLeft = highlight.offsetLeft
}

function get_selected_set() {
    const select = document.querySelector('select[name=re-set]')
    let set = select.options[select.selectedIndex].id
    if (set == 'selected-res') {
        return Array.from(document.querySelectorAll('#restriction-enzymes input:checked')).map(el => el.id)
    } else {
        return subsets[set]
    }
}
function show_re_set() {
    document.querySelectorAll('#restriction-enzymes tr').forEach(el => el.hidden = true)
    get_selected_set().forEach(name => document.getElementById(name).parentElement.parentElement.hidden = false)
    function filter_by(input, key) {
        const filter = new RegExp(input.value, 'i')
        document.querySelectorAll('#restriction-enzymes input').forEach(el => {
            if (!key(el).match(filter)) {
                el.parentElement.parentElement.hidden = true
            }
        })
    }
    filter_by(document.querySelector('#name-filter'), el => el.id)
    filter_by(document.querySelector('#seq-filter'), el => restriction_enzymes[el.id])
}
function select_shown_res(checked) {
    get_selected_set().forEach(name => document.getElementById(name).checked = checked)
}
function clear_re_set() {
    document.querySelectorAll('#restriction-enzymes input:checked').forEach(el => el.checked = false)
}
function reset_filters() {
    document.querySelectorAll('#name-filter, #seq-filter').forEach(el => {
        el.value = ''
        el.dispatchEvent(new Event('input'))
    })
}

let crosshair = document.getElementById('crosshair')
let crosshair_w = 0
crosshair.style.width = crosshair_w + 'px'

let drag_start = 0
let dragging = false

document.querySelectorAll('body, textarea').forEach(el => {
    el.addEventListener('keypress', e => {
        if (e.code == 'Enter' && e.ctrlKey) {
            e.preventDefault()
            calculate()
        }
    })
})
// document.querySelector('.re-select').addEventListener('click', e => {
//     if (e.button != 0) return;
//     document.getElementById('re-checkbox').checked = true
// })
document.querySelectorAll('.re-select').forEach(
    el => el.addEventListener('click', e => {
        if (e.button != 0) return;
        document.getElementById('re-checkbox').checked = true
    })
)
document.getElementById('results').addEventListener('mousemove', function (e) {
    if (dragging) {
        let left = Math.min(e.clientX, drag_start)
        let right = Math.max(e.clientX, drag_start)
        crosshair.style.left = left - crosshair_w / 2 + 'px'
        crosshair.style.width = right - left + crosshair_w + 'px'
    } else {
        crosshair.style.left = e.clientX - crosshair_w / 2 + 'px'
    }
    crosshair.style.top = document.getElementById('results').offsetTop + 'px'
    crosshair.style.height = document.getElementById('results').offsetHeight + 'px'
})
document.getElementById('results').addEventListener('mouseover', function (e) {
    crosshair.style.display = 'block'
})
document.getElementById('results').addEventListener('mouseout', function (e) {
    crosshair.style.display = 'none'
})
document.getElementById('results').addEventListener('mousedown', function (e) {
    if (e.button != 0) return;
    drag_start = e.clientX
    dragging = true
})
document.getElementById('results').addEventListener('mouseup', function (e) {
    if (e.button != 0) return;
    dragging = false
    crosshair.style.width = crosshair_w + 'px'
})