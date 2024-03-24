
function table_header(seq) {
    const clustalX_colorscheme = {
        '*': '#ffffff',
        'A': '#80a0f0', 'D': '#c048c0', 'W': '#80a0f0', 'C': '#f08080',
        'I': '#80a0f0', 'N': '#15c015', 'V': '#80a0f0', 'G': '#f09048',
        'L': '#80a0f0', 'Q': '#15c015', 'K': '#f01505', 'P': '#c0c000',
        'M': '#80a0f0', 'S': '#15c015', 'R': '#f01505', 'H': '#15a4a4',
        'F': '#80a0f0', 'T': '#15c015', 'E': '#c048c0', 'Y': '#15a4a4',
    }
    let numbers = ''
    let translated = translate(seq)
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

function result_tr(seq, re_name, fit) {
    function highlight_site(seq, pos, site) {
        return seq.substring(0, pos) +
            '<span class="highlight">' +
            site +
            '</span>' +
            seq.substring(pos + site.length)
    }
    let tr = document.createElement('tr')
    tr.innerHTML = 
        `<td><span class="action-link" onclick="scroll_to_match(this)" title="${restriction_enzymes[re_name]}">` + re_name + '</span></td>' +
        '<td>' + highlight_site(seq, fit.pos, fit.seq) + '</td>'
    return tr
}

function get_input_filtered(id, filter) {
    let input = document.getElementById('in-seq').value.toUpperCase()
    return input.split('').filter(c => filter.includes(c)).join('')
}
function get_input_cds(id) {
    return get_input_filtered(id, 'AGTC')
}
function get_input_prot(id) {
    return get_input_filtered(id, 'ACDEFGHIKLMNPQRSTVWY')
}

function calculate() {
    let input_type = document.querySelector('input[name="input-type"]:checked').id
    let cds = ''
    if (input_type == 'in-prot-radio') {
        const prot = get_input_prot('in-seq')
        cds = back_translate(prot)
    } else {
        cds = get_input_cds('in-seq')
        cds = cds.substring(0, cds.length - cds.length % 3)
    }

    let add_rev_comp = true
    let patterns = {}
    if (document.querySelector('#custom-seq-radio').checked) 
    {
        let patterns_in = document.querySelector('#custom-seq').value
            .toUpperCase().split('\n')
            .map(line => line.split('').filter(c => 'AGTCRYSWKMBVDHN'.includes(c)).join(''))
            .filter(pattern => pattern != '')
        patterns = Object.fromEntries(patterns_in.map(p => [p, p]))
        add_rev_comp = false
    } else {
        // get selected REs
        patterns = Object.fromEntries(Array.from(document.querySelectorAll('#restriction-enzymes input:checked')).map(el => [el.id, restriction_enzymes[el.id]]))
    }
    
    let results_table = document.getElementById('results')
    results_table.innerHTML = ''
    results_table.appendChild(table_header(cds))

    let found_seqs = []
    for (const pattern_name of Object.keys(patterns)) {
        fits = try_insert_pattern(cds, patterns[pattern_name], add_rev_comp)    
        if (fits.length > 0 && !found_seqs.includes(patterns[pattern_name])) {
            //add finds to table
            fits.map(fit => result_tr(cds, pattern_name, fit)).forEach(tr => document.querySelector('table').appendChild(tr));
            found_seqs.push(patterns[pattern_name])
        }
    }
    if (found_seqs.length == 0) {
        results_table.innerHTML += '<p style="color: red; position: sticky; left: 0;">Nothing found :(</p>'
    }
    results_table.scrollLeft = 0
    results_table.scrollTop = 0
}

function rev_comp_input() {
    document.getElementById('in-seq').value = rev_comp(get_input_cds())
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
function select_subset(el) {
    subsets[el.dataset.subset].forEach(
        name => document.querySelector('#restriction-enzymes #' + name).checked = true
    )
}
function deselect_subset(el) {
    subsets[el.dataset.subset].forEach(name => document.querySelector('#restriction-enzymes #' + name).checked = false)
}

let crosshair = document.getElementById('crosshair')
let crosshair_w = 0
crosshair.style.width = crosshair_w + 'px'

let drag_start = 0
let dragging = false

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