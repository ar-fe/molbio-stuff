function rev_comp(seq) {
    const comp = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'
    }
    return seq.split('').reverse().map(c => comp[c]).join('')
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

let REs_raw = {
    aani: {
        name: "AanI",
        seq: "TTATAA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    aari: {
        name: "AarI",
        seq: "CACCTGC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    aatii: {
        name: "AatII",
        seq: "GACGTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    acc65i: {
        name: "Acc65I",
        seq: "GGTACC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    adei: {
        name: "AdeI",
        seq: "CACNNNGTG",
        cuts: [6],
        dam: false,
        dcm: false
    },
    ajii: {
        name: "AjiI",
        seq: "CACGTC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    ajui: {
        name: "AjuI",
        seq: "GAA(N)7TTGG",
        cuts: [-7,
            11],
        dam: false,
        dcm: false
    },
    aloi: {
        name: "AloI",
        seq: "GAAC(N)6TCC",
        cuts: [-7,
            12],
        dam: false,
        dcm: false
    },
    alui: {
        name: "AluI",
        seq: "AGCT",
        cuts: [2],
        dam: false,
        dcm: false
    },
    alw21i: {
        name: "Alw21I",
        seq: "GWGCWC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    alw26i: {
        name: "Alw26I",
        seq: "GTCTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    alw44i: {
        name: "Alw44I",
        seq: "GTGCAC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    apai: {
        name: "ApaI",
        seq: "GGGCCC",
        cuts: [5],
        dam: false,
        dcm: true
    },
    bamhi: {
        name: "BamHI",
        seq: "GGATCC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    baui: {
        name: "BauI",
        seq: "CACGAG",
        cuts: [-5],
        dam: false,
        dcm: false
    },
    bcli: {
        name: "BclI",
        seq: "TGATCA",
        cuts: [1],
        dam: true,
        dcm: false
    },
    bcni: {
        name: "BcnI",
        seq: "CCSGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bcui: {
        name: "BcuI",
        seq: "ACTAGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bfmi: {
        name: "BfmI",
        seq: "CTRYAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bfui: {
        name: "BfuI",
        seq: "GTATCC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    bgli: {
        name: "BglI",
        seq: "GCCNNNNNGGC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    bglii: {
        name: "BglII",
        seq: "AGATCT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bme1390i: {
        name: "Bme1390I",
        seq: "CCNGG",
        cuts: [2],
        dam: false,
        dcm: true
    },
    boxi: {
        name: "BoxI",
        seq: "GACNNNNGTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bpii: {
        name: "BpiI",
        seq: "GAAGAC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bpli: {
        name: "BplI",
        seq: "GAG(N)5CTC",
        cuts: [-8,
            13],
        dam: false,
        dcm: false
    },
    bpu10i: {
        name: "Bpu10I",
        seq: "CCTNAGC",
        cuts: [-5],
        dam: false,
        dcm: false
    },
    bpu1102i: {
        name: "Bpu1102I",
        seq: "GCTNAGC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsedi: {
        name: "BseDI",
        seq: "CCNNGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsegi: {
        name: "BseGI",
        seq: "GGATG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bseji: {
        name: "BseJI",
        seq: "GATNNNNATC",
        cuts: [5],
        dam: true,
        dcm: false
    },
    bseli: {
        name: "BseLI",
        seq: "CCNNNNNNNGG",
        cuts: [7],
        dam: false,
        dcm: true
    },
    bsemi: {
        name: "BseMI",
        seq: "GCAATG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsemii: {
        name: "BseMII",
        seq: "CTCAG",
        cuts: [10],
        dam: false,
        dcm: false
    },
    bseni: {
        name: "BseNI",
        seq: "ACTGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsesi: {
        name: "BseSI",
        seq: "GKGCMC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bsexi: {
        name: "BseXI",
        seq: "GCAGC",
        cuts: [8],
        dam: false,
        dcm: false
    },
    bsh1236i: {
        name: "Bsh1236I",
        seq: "CGCG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsh1285i: {
        name: "Bsh1285I",
        seq: "CGRYCG",
        cuts: [4],
        dam: true,
        dcm: false
    },
    bshni: {
        name: "BshNI",
        seq: "GGYRCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    bshti: {
        name: "BshTI",
        seq: "ACCGGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsp119i: {
        name: "Bsp119I",
        seq: "TTCGAA",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsp120i: {
        name: "Bsp120I",
        seq: "GGGCCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    bsp1407i: {
        name: "Bsp1407I",
        seq: "TGTACA",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsp143i: {
        name: "Bsp143I",
        seq: "GATC",
        cuts: [0],
        dam: false,
        dcm: false
    },
    bsp68i: {
        name: "Bsp68I",
        seq: "TCGCGA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    bspli: {
        name: "BspLI",
        seq: "GGNNCC",
        cuts: [3],
        dam: false,
        dcm: true
    },
    bspoi: {
        name: "BspOI",
        seq: "GCTAGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bsppi: {
        name: "BspPI",
        seq: "GGATC",
        cuts: [4],
        dam: true,
        dcm: false
    },
    bspti: {
        name: "BspTI",
        seq: "CTTAAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bst1107i: {
        name: "Bst1107I",
        seq: "GTATAC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    bstxi: {
        name: "BstXI",
        seq: "CCANNNNNNTGG",
        cuts: [8],
        dam: false,
        dcm: true
    },
    bsu15i: {
        name: "Bsu15I",
        seq: "ATCGAT",
        cuts: [2],
        dam: true,
        dcm: false
    },
    bsuri: {
        name: "BsuRI",
        seq: "GGCC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bvei: {
        name: "BveI",
        seq: "ACCTGC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    caii: {
        name: "CaiI",
        seq: "CAGNNNCTG",
        cuts: [6],
        dam: false,
        dcm: true
    },
    cfr10i: {
        name: "Cfr10I",
        seq: "RCCGGY",
        cuts: [1],
        dam: false,
        dcm: false
    },
    cfr13i: {
        name: "Cfr13I",
        seq: "GGNCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    cfr42i: {
        name: "Cfr42I",
        seq: "CCGCGG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    cfr9i: {
        name: "Cfr9I",
        seq: "CCCGGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    cpoi: {
        name: "CpoI",
        seq: "CGGWCCG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    csei: {
        name: "CseI",
        seq: "GACGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    csp6i: {
        name: "Csp6I",
        seq: "GTAC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    dpni: {
        name: "DpnI",
        seq: "GATC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    drai: {
        name: "DraI",
        seq: "TTTAAA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    eam1104i: {
        name: "Eam1104I",
        seq: "CTCTTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    eam1105i: {
        name: "Eam1105I",
        seq: "GACNNNNNGTC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    ecl136ii: {
        name: "Ecl136II",
        seq: "GAGCTC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    eco105i: {
        name: "Eco105I",
        seq: "TACGTA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    eco130i: {
        name: "Eco130I",
        seq: "CCWWGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    eco147i: {
        name: "Eco147I",
        seq: "AGGCCT",
        cuts: [3],
        dam: false,
        dcm: true
    },
    eco24i: {
        name: "Eco24I",
        seq: "GRGCYC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    eco31i: {
        name: "Eco31I",
        seq: "GGTCTC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    eco32i: {
        name: "Eco32I",
        seq: "GATATC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    eco47i: {
        name: "Eco47I",
        seq: "GGWCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    eco47iii: {
        name: "Eco47III",
        seq: "AGCGCT",
        cuts: [3],
        dam: false,
        dcm: false
    },
    eco52i: {
        name: "Eco52I",
        seq: "CGGCCG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    eco57i: {
        name: "Eco57I",
        seq: "CTGAAG",
        cuts: [16],
        dam: false,
        dcm: false
    },
    eco72i: {
        name: "Eco72I",
        seq: "CACGTG",
        cuts: [3],
        dam: false,
        dcm: false
    },
    eco81i: {
        name: "Eco81I",
        seq: "CCTNAGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    eco88i: {
        name: "Eco88I",
        seq: "CYCGRG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    eco91i: {
        name: "Eco91I",
        seq: "GGTNACC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    ecoo109i: {
        name: "EcoO109I",
        seq: "RGGNCCY",
        cuts: [2],
        dam: false,
        dcm: true
    },
    ecori: {
        name: "EcoRI",
        seq: "GAATTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    ecorii: {
        name: "EcoRII",
        seq: "CCWGG",
        cuts: [0],
        dam: false,
        dcm: true
    },
    ehei: {
        name: "EheI",
        seq: "GGCGCC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    esp3i: {
        name: "Esp3I",
        seq: "CGTCTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    faqi: {
        name: "FaqI",
        seq: "GGGAC",
        cuts: [10],
        dam: false,
        dcm: false
    },
    fspai: {
        name: "FspAI",
        seq: "RTGCGCAY",
        cuts: [4],
        dam: false,
        dcm: false
    },
    fspbi: {
        name: "FspBI",
        seq: "CTAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    gsui: {
        name: "GsuI",
        seq: "CTGGAG",
        cuts: [16],
        dam: false,
        dcm: true
    },
    hhai: {
        name: "HhaI",
        seq: "GCGC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    hin1i: {
        name: "Hin1I",
        seq: "GRCGYC",
        cuts: [2],
        dam: false,
        dcm: true
    },
    hin1ii: {
        name: "Hin1II",
        seq: "CATG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    hin6i: {
        name: "Hin6I",
        seq: "GCGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    hincii: {
        name: "HincII",
        seq: "GTYRAC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    hindiii: {
        name: "HindIII",
        seq: "AAGCTT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    hinfi: {
        name: "HinfI",
        seq: "GANTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    hpaii: {
        name: "HpaII",
        seq: "CCGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    hphi: {
        name: "HphI",
        seq: "GGTGA",
        cuts: [8],
        dam: true,
        dcm: true
    },
    hpy8i: {
        name: "Hpy8I",
        seq: "GTNNAC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    hpyf10vi: {
        name: "HpyF10VI",
        seq: "GCNNNNNNNGC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    hpyf3i: {
        name: "HpyF3I",
        seq: "CTNAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    kpn2i: {
        name: "Kpn2I",
        seq: "TCCGGA",
        cuts: [1],
        dam: false,
        dcm: false
    },
    kpni: {
        name: "KpnI",
        seq: "GGTACC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    kspai: {
        name: "KspAI",
        seq: "GTTAAC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    lgui: {
        name: "LguI",
        seq: "GCTCTTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    lsp1109i: {
        name: "Lsp1109I",
        seq: "GCAGC",
        cuts: [8],
        dam: false,
        dcm: false
    },
    lwei: {
        name: "LweI",
        seq: "GCATC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    maubi: {
        name: "MauBI",
        seq: "CGCGCGCG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    mbii: {
        name: "MbiI",
        seq: "CCGCTC",
        cuts: [-3],
        dam: false,
        dcm: false
    },
    mboi: {
        name: "MboI",
        seq: "GATC",
        cuts: [0],
        dam: true,
        dcm: false
    },
    mboii: {
        name: "MboII",
        seq: "GAAGA",
        cuts: [8],
        dam: true,
        dcm: false
    },
    mlsi: {
        name: "MlsI",
        seq: "TGGCCA",
        cuts: [3],
        dam: false,
        dcm: true
    },
    mlui: {
        name: "MluI",
        seq: "ACGCGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    mnli: {
        name: "MnlI",
        seq: "CCTC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    mph1103i: {
        name: "Mph1103I",
        seq: "ATGCAT",
        cuts: [5],
        dam: false,
        dcm: false
    },
    mrei: {
        name: "MreI",
        seq: "CGCCGGCG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    mspi: {
        name: "MspI",
        seq: "CCGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    mssi: {
        name: "MssI",
        seq: "GTTTAAAC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    muni: {
        name: "MunI",
        seq: "CAATTG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    mva1269i: {
        name: "Mva1269I",
        seq: "GAATGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    mvai: {
        name: "MvaI",
        seq: "CCWGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    ncoi: {
        name: "NcoI",
        seq: "CCATGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    ndei: {
        name: "NdeI",
        seq: "CATATG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    nhei: {
        name: "NheI",
        seq: "GCTAGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    nmuci: {
        name: "NmuCI",
        seq: "GTSAC",
        cuts: [0],
        dam: false,
        dcm: false
    },
    noti: {
        name: "NotI",
        seq: "GCGGCCGC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    nsbi: {
        name: "NsbI",
        seq: "TGCGCA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    olii: {
        name: "OliI",
        seq: "CACNNNNGTG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    paci: {
        name: "PacI",
        seq: "TTAATTAA",
        cuts: [5],
        dam: false,
        dcm: false
    },
    paei: {
        name: "PaeI",
        seq: "GCATGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    pagi: {
        name: "PagI",
        seq: "TCATGA",
        cuts: [1],
        dam: true,
        dcm: false
    },
    pasi: {
        name: "PasI",
        seq: "CCCWGGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    paui: {
        name: "PauI",
        seq: "GCGCGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    pdii: {
        name: "PdiI",
        seq: "GCCGGC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    pdmi: {
        name: "PdmI",
        seq: "GAANNNNTTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    pfei: {
        name: "PfeI",
        seq: "GAWTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    pfl23ii: {
        name: "Pfl23II",
        seq: "CGTACG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    pfoi: {
        name: "PfoI",
        seq: "TCCNGGA",
        cuts: [1],
        dam: true,
        dcm: true
    },
    ppu21i: {
        name: "Ppu21I",
        seq: "YACGTR",
        cuts: [3],
        dam: false,
        dcm: false
    },
    psci: {
        name: "PscI",
        seq: "ACATGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    psp1406i: {
        name: "Psp1406I",
        seq: "AACGTT",
        cuts: [2],
        dam: false,
        dcm: false
    },
    psp5ii: {
        name: "Psp5II",
        seq: "RGGWCCY",
        cuts: [2],
        dam: false,
        dcm: true
    },
    psti: {
        name: "PstI",
        seq: "CTGCAG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    psui: {
        name: "PsuI",
        seq: "RGATCY",
        cuts: [1],
        dam: false,
        dcm: false
    },
    psyi: {
        name: "PsyI",
        seq: "GACNNNGTC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    pvui: {
        name: "PvuI",
        seq: "CGATCG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    pvuii: {
        name: "PvuII",
        seq: "CAGCTG",
        cuts: [3],
        dam: false,
        dcm: false
    },
    rsai: {
        name: "RsaI",
        seq: "GTAC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    rsei: {
        name: "RseI",
        seq: "CAYNNNNRTG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    saci: {
        name: "SacI",
        seq: "GAGCTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    sali: {
        name: "SalI",
        seq: "GTCGAC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    sati: {
        name: "SatI",
        seq: "GCNGC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    scai: {
        name: "ScaI",
        seq: "AGTACT",
        cuts: [3],
        dam: false,
        dcm: false
    },
    schi: {
        name: "SchI",
        seq: "GAGTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    sdai: {
        name: "SdaI",
        seq: "CCTGCAGG",
        cuts: [6],
        dam: false,
        dcm: false
    },
    sdui: {
        name: "SduI",
        seq: "GDGCHC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    sfaai: {
        name: "SfaAI",
        seq: "GCGATCGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    sfii: {
        name: "SfiI",
        seq: "GGCCNNNNNGGCC",
        cuts: [8],
        dam: false,
        dcm: true
    },
    sgei: {
        name: "SgeI",
        seq: "m5CNNG",
        cuts: [9],
        dam: false,
        dcm: false
    },
    sgrdi: {
        name: "SgrDI",
        seq: "CGTCGACG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    sgsi: {
        name: "SgsI",
        seq: "GGCGCGCC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    smai: {
        name: "SmaI",
        seq: "CCCGGG",
        cuts: [3],
        dam: false,
        dcm: false
    },
    smii: {
        name: "SmiI",
        seq: "ATTTAAAT",
        cuts: [4],
        dam: false,
        dcm: false
    },
    smoi: {
        name: "SmoI",
        seq: "CTYRAG",
        cuts: [1],
        dam: true,
        dcm: false
    },
    ssii: {
        name: "SsiI",
        seq: "CCGC",
        cuts: [-3],
        dam: false,
        dcm: false
    },
    sspdi: {
        name: "SspDI",
        seq: "GGCGCC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    sspi: {
        name: "SspI",
        seq: "AATATT",
        cuts: [3],
        dam: false,
        dcm: false
    },
    taai: {
        name: "TaaI",
        seq: "ACNGT",
        cuts: [3],
        dam: false,
        dcm: false
    },
    taii: {
        name: "TaiI",
        seq: "ACGT",
        cuts: [4],
        dam: false,
        dcm: false
    },
    taqi: {
        name: "TaqI",
        seq: "TCGA",
        cuts: [1],
        dam: true,
        dcm: false
    },
    tasi: {
        name: "TasI",
        seq: "AATT",
        cuts: [0],
        dam: false,
        dcm: false
    },
    tati: {
        name: "TatI",
        seq: "WGTACW",
        cuts: [1],
        dam: false,
        dcm: false
    },
    taui: {
        name: "TauI",
        seq: "GCSGC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    tru1i: {
        name: "Tru1I",
        seq: "TTAA",
        cuts: [1],
        dam: false,
        dcm: false
    },
    tscai: {
        name: "TscAI",
        seq: "CASTG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    van91i: {
        name: "Van91I",
        seq: "CCANNNNNTGG",
        cuts: [7],
        dam: false,
        dcm: true
    },
    vspi: {
        name: "VspI",
        seq: "ATTAAT",
        cuts: [2],
        dam: false,
        dcm: false
    },
    xagi: {
        name: "XagI",
        seq: "CCTNNNNNAGG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    xapi: {
        name: "XapI",
        seq: "RAATTY",
        cuts: [1],
        dam: false,
        dcm: false
    },
    xbai: {
        name: "XbaI",
        seq: "TCTAGA",
        cuts: [1],
        dam: true,
        dcm: false
    },
    xcei: {
        name: "XceI",
        seq: "RCATGY",
        cuts: [5],
        dam: false,
        dcm: false
    },
    xhoi: {
        name: "XhoI",
        seq: "CTCGAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    xmaji: {
        name: "XmaJI",
        seq: "CCTAGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    xmii: {
        name: "XmiI",
        seq: "GTMKAC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    absi: {
        name: "AbsI",
        seq: "CCTCGAGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    acc16i: {
        name: "Acc16I",
        seq: "TGCGCA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    acc36i: {
        name: "Acc36I",
        seq: "ACCTGC",
        cuts: [10],
        dam: false,
        dcm: false
    },
    accb1i: {
        name: "AccB1I",
        seq: "GGYRCC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    accb7i: {
        name: "AccB7I",
        seq: "CCANNNNNTGG",
        cuts: [7],
        dam: false,
        dcm: true
    },
    accbsi: {
        name: "AccBSI",
        seq: "GAGCGG",
        cuts: [3],
        dam: false,
        dcm: false
    },
    acli: {
        name: "AclI",
        seq: "AACGTT",
        cuts: [2],
        dam: false,
        dcm: false
    },
    aclwi: {
        name: "AclWI",
        seq: "GGATC",
        cuts: [9],
        dam: true,
        dcm: false
    },
    acoi: {
        name: "AcoI",
        seq: "YGGCCR",
        cuts: [1],
        dam: false,
        dcm: true
    },
    acsi: {
        name: "AcsI",
        seq: "RAATTY",
        cuts: [1],
        dam: false,
        dcm: false
    },
    acui: {
        name: "AcuI",
        seq: "CTGAAG",
        cuts: [16],
        dam: false,
        dcm: false
    },
    afei: {
        name: "AfeI",
        seq: "AGCGCT",
        cuts: [3],
        dam: false,
        dcm: false
    },
    agsi: {
        name: "AgsI",
        seq: "TTSAA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    ahli: {
        name: "AhlI",
        seq: "ACTAGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    ajni: {
        name: "AjnI",
        seq: "CCWGG",
        cuts: [0],
        dam: false,
        dcm: false
    },
    alei: {
        name: "AleI",
        seq: "CACNNNNGTG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    alubi: {
        name: "AluBI",
        seq: "AGCT",
        cuts: [2],
        dam: false,
        dcm: false
    },
    ama87i: {
        name: "Ama87I",
        seq: "CYCGRG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    arsi: {
        name: "ArsI",
        seq: "GACNNNNNNTTYG",
        cuts: [-8,
            24],
        dam: false,
        dcm: false
    },
    asigi: {
        name: "AsiGI",
        seq: "ACCGGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    asisi: {
        name: "AsiSI",
        seq: "GCGATCGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    aspa2i: {
        name: "AspA2I",
        seq: "CCTAGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    asplei: {
        name: "AspLEI",
        seq: "GCGC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    asps9i: {
        name: "AspS9I",
        seq: "GGNCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    asuc2i: {
        name: "AsuC2I",
        seq: "CCSGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    asuhpi: {
        name: "AsuHPI",
        seq: "GGTGA",
        cuts: [13],
        dam: true,
        dcm: false
    },
    asunhi: {
        name: "AsuNHI",
        seq: "GCTAGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bari: {
        name: "BarI",
        seq: "GAAGNNNNNNTAC",
        cuts: [-7,
            25],
        dam: false,
        dcm: false
    },
    bbv12i: {
        name: "Bbv12I",
        seq: "GWGCWC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bme18i: {
        name: "Bme18I",
        seq: "GGWCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    bmoi: {
        name: "BmoI",
        seq: "GGGWCCC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bmti: {
        name: "BmtI",
        seq: "GCTAGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bmui: {
        name: "BmuI",
        seq: "ACTGGG",
        cuts: [11],
        dam: false,
        dcm: false
    },
    bpmi: {
        name: "BpmI",
        seq: "CTGGAG",
        cuts: [16],
        dam: false,
        dcm: false
    },
    bpu14i: {
        name: "Bpu14I",
        seq: "TTCGAA",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsa29i: {
        name: "Bsa29I",
        seq: "ATCGAT",
        cuts: [2],
        dam: true,
        dcm: false
    },
    bsc4i: {
        name: "Bsc4I",
        seq: "CCNNNNNNNGG",
        cuts: [7],
        dam: false,
        dcm: false
    },
    bse1i: {
        name: "Bse1I",
        seq: "ACTGG",
        cuts: [6],
        dam: false,
        dcm: false
    },
    bse118i: {
        name: "Bse118I",
        seq: "RCCGGY",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bse21i: {
        name: "Bse21I",
        seq: "CCTNAGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bse3di: {
        name: "Bse3DI",
        seq: "GCAATG",
        cuts: [8],
        dam: false,
        dcm: false
    },
    bse8i: {
        name: "Bse8I",
        seq: "GATNNNNATC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bsepi: {
        name: "BsePI",
        seq: "GCGCGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsex3i: {
        name: "BseX3I",
        seq: "CGGCCG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bslfi: {
        name: "BslFI",
        seq: "GGGAC",
        cuts: [15],
        dam: false,
        dcm: false
    },
    bso31i: {
        name: "Bso31I",
        seq: "GGTCTC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    bsp13i: {
        name: "Bsp13I",
        seq: "TCCGGA",
        cuts: [1],
        dam: true,
        dcm: false
    },
    bsp1720i: {
        name: "Bsp1720I",
        seq: "GCTNAGC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsp19i: {
        name: "Bsp19I",
        seq: "CCATGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bspaci: {
        name: "BspACI",
        seq: "CCGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bspfni: {
        name: "BspFNI",
        seq: "CGCG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsseci: {
        name: "BssECI",
        seq: "CCNNGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bssnai: {
        name: "BssNAI",
        seq: "GTATAC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    bsst1i: {
        name: "BssT1I",
        seq: "CCWWGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bst2bi: {
        name: "Bst2BI",
        seq: "CTCGTG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bst2ui: {
        name: "Bst2UI",
        seq: "CCWGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bst4ci: {
        name: "Bst4CI",
        seq: "ACNGT",
        cuts: [3],
        dam: false,
        dcm: false
    },
    bst6i: {
        name: "Bst6I",
        seq: "CTCTTC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    bstaci: {
        name: "BstACI",
        seq: "GRCGYC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bstafi: {
        name: "BstAFI",
        seq: "CTTAAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bstapi: {
        name: "BstAPI",
        seq: "GCANNNNNTGC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    bstaui: {
        name: "BstAUI",
        seq: "TGTACA",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bstbai: {
        name: "BstBAI",
        seq: "YACGTR",
        cuts: [3],
        dam: false,
        dcm: false
    },
    bstc8i: {
        name: "BstC8I",
        seq: "GCNNGC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    bstdei: {
        name: "BstDEI",
        seq: "CTNAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bstdsi: {
        name: "BstDSI",
        seq: "CCRYGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsteni: {
        name: "BstENI",
        seq: "CCTNNNNNAGG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bstf5i: {
        name: "BstF5I",
        seq: "GGATG",
        cuts: [7],
        dam: false,
        dcm: false
    },
    bstfni: {
        name: "BstFNI",
        seq: "CGCG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsth2i: {
        name: "BstH2I",
        seq: "RGCGCY",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bsthhi: {
        name: "BstHHI",
        seq: "GCGC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    bstkti: {
        name: "BstKTI",
        seq: "GATC",
        cuts: [3],
        dam: true,
        dcm: false
    },
    bstmai: {
        name: "BstMAI",
        seq: "GTCTC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    bstmbi: {
        name: "BstMBI",
        seq: "GATC",
        cuts: [0],
        dam: true,
        dcm: false
    },
    bstmci: {
        name: "BstMCI",
        seq: "CGRYCG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    bstmwi: {
        name: "BstMWI",
        seq: "GCNNNNNNNGC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    bstnsi: {
        name: "BstNSI",
        seq: "RCATGY",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bstpai: {
        name: "BstPAI",
        seq: "GACNNNNGTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bstsci: {
        name: "BstSCI",
        seq: "CCNGG",
        cuts: [0],
        dam: false,
        dcm: true
    },
    bstsfi: {
        name: "BstSFI",
        seq: "CTRYAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bstsli: {
        name: "BstSLI",
        seq: "GKGCMC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bstsni: {
        name: "BstSNI",
        seq: "TACGTA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    bstv1i: {
        name: "BstV1I",
        seq: "GCAGC",
        cuts: [13],
        dam: false,
        dcm: false
    },
    bstv2i: {
        name: "BstV2I",
        seq: "GAAGAC",
        cuts: [8],
        dam: false,
        dcm: false
    },
    bstx2i: {
        name: "BstX2I",
        seq: "RGATCY",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsui: {
        name: "BsuI",
        seq: "GTATCC",
        cuts: [12],
        dam: false,
        dcm: false
    },
    btri: {
        name: "BtrI",
        seq: "CACGTC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    ccii: {
        name: "CciI",
        seq: "TCATGA",
        cuts: [1],
        dam: false,
        dcm: false
    },
    ccini: {
        name: "CciNI",
        seq: "GCGGCCGC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    draiii: {
        name: "DraIII",
        seq: "CACNNNGTG",
        cuts: [6],
        dam: false,
        dcm: false
    },
    drii: {
        name: "DriI",
        seq: "GACNNNNNGTC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    dsedi: {
        name: "DseDI",
        seq: "GACNNNNNNGTC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    ecoicri: {
        name: "EcoICRI",
        seq: "GAGCTC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    ecorv: {
        name: "EcoRV",
        seq: "GATATC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    egei: {
        name: "EgeI",
        seq: "GGCGCC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    erhi: {
        name: "ErhI",
        seq: "CCWWGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    faei: {
        name: "FaeI",
        seq: "CATG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    faii: {
        name: "FaiI",
        seq: "YATR",
        cuts: [2],
        dam: false,
        dcm: false
    },
    fali: {
        name: "FalI",
        seq: "AAGNNNNNCTT",
        cuts: [-8,
            24],
        dam: false,
        dcm: false
    },
    fati: {
        name: "FatI",
        seq: "CATG",
        cuts: [0],
        dam: false,
        dcm: false
    },
    faui: {
        name: "FauI",
        seq: "CCCGC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    faundi: {
        name: "FauNDI",
        seq: "CATATG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    fbli: {
        name: "FblI",
        seq: "GTMKAC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    foki: {
        name: "FokI",
        seq: "GGATG",
        cuts: [9],
        dam: false,
        dcm: true
    },
    frioi: {
        name: "FriOI",
        seq: "GRGCYC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    fsp4hi: {
        name: "Fsp4HI",
        seq: "GCNGC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    gsai: {
        name: "GsaI",
        seq: "CCCAGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    haeiii: {
        name: "HaeIII",
        seq: "GGCC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    hgai: {
        name: "HgaI",
        seq: "GACGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    hindii: {
        name: "HindII",
        seq: "GTYRAC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    hpai: {
        name: "HpaI",
        seq: "GTTAAC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    hpyse526i: {
        name: "HpySE526I",
        seq: "ACGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    hspai: {
        name: "HspAI",
        seq: "GCGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    kroni: {
        name: "KroNI",
        seq: "GCCGGC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    ksp22i: {
        name: "Ksp22I",
        seq: "TGATCA",
        cuts: [1],
        dam: true,
        dcm: false
    },
    kzo9i: {
        name: "Kzo9I",
        seq: "GATC",
        cuts: [0],
        dam: false,
        dcm: false
    },
    lmni: {
        name: "LmnI",
        seq: "GCTCC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    mabi: {
        name: "MabI",
        seq: "ACCWGGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    mfei: {
        name: "MfeI",
        seq: "CAATTG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    mhli: {
        name: "MhlI",
        seq: "GDGCHC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    mly113i: {
        name: "Mly113I",
        seq: "GGCGCC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    mox20i: {
        name: "Mox20I",
        seq: "TGGCCA",
        cuts: [3],
        dam: false,
        dcm: true
    },
    mroni: {
        name: "MroNI",
        seq: "GCCGGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    mroxi: {
        name: "MroXI",
        seq: "GAANNNNTTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    mspa1i: {
        name: "MspA1I",
        seq: "CMGCKG",
        cuts: [3],
        dam: false,
        dcm: false
    },
    mspr9i: {
        name: "MspR9I",
        seq: "CCNGG",
        cuts: [2],
        dam: false,
        dcm: true
    },
    nrui: {
        name: "NruI",
        seq: "TCGCGA",
        cuts: [3],
        dam: true,
        dcm: false
    },
    palai: {
        name: "PalAI",
        seq: "GGCGCGCC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    pcei: {
        name: "PceI",
        seq: "AGGCCT",
        cuts: [3],
        dam: false,
        dcm: false
    },
    pcii: {
        name: "PciI",
        seq: "ACATGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    pcisi: {
        name: "PciSI",
        seq: "GCTCTTC",
        cuts: [8],
        dam: false,
        dcm: false
    },
    pcti: {
        name: "PctI",
        seq: "GAATGC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    peci: {
        name: "PecI",
        seq: "TTATAA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    ple19i: {
        name: "Ple19I",
        seq: "CGATCG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    ppsi: {
        name: "PpsI",
        seq: "GAGTC",
        cuts: [9],
        dam: false,
        dcm: false
    },
    pse31i: {
        name: "Pse31I",
        seq: "GGTCTC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    psii: {
        name: "PsiI",
        seq: "TTATAA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    psp124bi: {
        name: "Psp124BI",
        seq: "GAGCTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    psp6i: {
        name: "Psp6I",
        seq: "CCWGG",
        cuts: [0],
        dam: true,
        dcm: false
    },
    pspci: {
        name: "PspCI",
        seq: "CACGTG",
        cuts: [3],
        dam: false,
        dcm: false
    },
    pspei: {
        name: "PspEI",
        seq: "GGTNACC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    pspli: {
        name: "PspLI",
        seq: "CGTACG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    pspn4i: {
        name: "PspN4I",
        seq: "GGNNCC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    pspomi: {
        name: "PspOMI",
        seq: "GGGCCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    pspppi: {
        name: "PspPPI",
        seq: "RGGWCCY",
        cuts: [2],
        dam: false,
        dcm: true
    },
    pspxi: {
        name: "PspXI",
        seq: "VCTCGAGB",
        cuts: [2],
        dam: false,
        dcm: false
    },
    psri: {
        name: "PsrI",
        seq: "GAACNNNNNNTAC",
        cuts: [-7,
            25],
        dam: false,
        dcm: false
    },
    pstni: {
        name: "PstNI",
        seq: "CAGNNNCTG",
        cuts: [6],
        dam: false,
        dcm: false
    },
    rgai: {
        name: "RgaI",
        seq: "GCGATCGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    rigi: {
        name: "RigI",
        seq: "GGCCGGCC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    rsani: {
        name: "RsaNI",
        seq: "GTAC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    rsr2i: {
        name: "Rsr2I",
        seq: "CGGWCCG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    sbfi: {
        name: "SbfI",
        seq: "CCTGCAGG",
        cuts: [6],
        dam: false,
        dcm: false
    },
    seti: {
        name: "SetI",
        seq: "ASST",
        cuts: [4],
        dam: false,
        dcm: false
    },
    sfani: {
        name: "SfaNI",
        seq: "GCATC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    sfr274i: {
        name: "Sfr274I",
        seq: "CTCGAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    sfr303i: {
        name: "Sfr303I",
        seq: "CCGCGG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    smimi: {
        name: "SmiMI",
        seq: "CAYNNNNRTG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    sphi: {
        name: "SphI",
        seq: "GCATGC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    sse9i: {
        name: "Sse9I",
        seq: "AATT",
        cuts: [0],
        dam: false,
        dcm: false
    },
    sspmi: {
        name: "SspMI",
        seq: "CTAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    tru9i: {
        name: "Tru9I",
        seq: "TTAA",
        cuts: [1],
        dam: false,
        dcm: false
    },
    tsefi: {
        name: "TseFI",
        seq: "GTSAC",
        cuts: [0],
        dam: false,
        dcm: false
    },
    tth111i: {
        name: "Tth111I",
        seq: "GACNNNGTC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    vnei: {
        name: "VneI",
        seq: "GTGCAC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    xmai: {
        name: "XmaI",
        seq: "CCCGGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    zrai: {
        name: "ZraI",
        seq: "GACGTC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    zrmi: {
        name: "ZrmI",
        seq: "AGTACT",
        cuts: [3],
        dam: false,
        dcm: false
    },
    zsp2i: {
        name: "Zsp2I",
        seq: "ATGCAT",
        cuts: [5],
        dam: false,
        dcm: false
    },
    acci: {
        name: "AccI",
        seq: "GTMKAC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    acii: {
        name: "AciI",
        seq: "CCGC",
        cuts: [-3],
        dam: false,
        dcm: false
    },
    aflii: {
        name: "AflII",
        seq: "CTTAAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    afliii: {
        name: "AflIII",
        seq: "ACRYGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    agei: {
        name: "AgeI",
        seq: "ACCGGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    ahdi: {
        name: "AhdI",
        seq: "GACNNNNNGTC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    alwi: {
        name: "AlwI",
        seq: "GGATC",
        cuts: [4],
        dam: true,
        dcm: false
    },
    alwni: {
        name: "AlwNI",
        seq: "CAGNNNCTG",
        cuts: [6],
        dam: false,
        dcm: true
    },
    apali: {
        name: "ApaLI",
        seq: "GTGCAC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    apeki: {
        name: "ApeKI",
        seq: "GCWGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    apoi: {
        name: "ApoI",
        seq: "RAATTY",
        cuts: [1],
        dam: false,
        dcm: false
    },
    asci: {
        name: "AscI",
        seq: "GGCGCGCC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    asei: {
        name: "AseI",
        seq: "ATTAAT",
        cuts: [2],
        dam: false,
        dcm: false
    },
    avai: {
        name: "AvaI",
        seq: "CYCGRG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    avaii: {
        name: "AvaII",
        seq: "GGWCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    avrii: {
        name: "AvrII",
        seq: "CCTAGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    baegi: {
        name: "BaeGI",
        seq: "GKGCMC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    baei: {
        name: "BaeI",
        seq: "",
        cuts: [10],
        dam: false,
        dcm: false
    },
    bani: {
        name: "BanI",
        seq: "GGYRCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    banii: {
        name: "BanII",
        seq: "GRGCYC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bbsi: {
        name: "BbsI",
        seq: "GAAGAC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bbvci: {
        name: "BbvCI",
        seq: "CCTCAGC",
        cuts: [-5],
        dam: false,
        dcm: false
    },
    bbvi: {
        name: "BbvI",
        seq: "GCAGC",
        cuts: [8],
        dam: false,
        dcm: false
    },
    bcci: {
        name: "BccI",
        seq: "CCATC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    bceai: {
        name: "BceAI",
        seq: "ACGGC",
        cuts: [12],
        dam: false,
        dcm: false
    },
    bcgi: {
        name: "BcgI",
        seq: "",
        cuts: [10],
        dam: true,
        dcm: false
    },
    bcivi: {
        name: "BciVI",
        seq: "GTATCC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    bcodi: {
        name: "BcoDI",
        seq: "GTCTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bfai: {
        name: "BfaI",
        seq: "CTAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bfuai: {
        name: "BfuAI",
        seq: "ACCTGC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    blpi: {
        name: "BlpI",
        seq: "GCTNAGC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bmgbi: {
        name: "BmgBI",
        seq: "CACGTC",
        cuts: [-3],
        dam: false,
        dcm: false
    },
    bmri: {
        name: "BmrI",
        seq: "ACTGGG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bpuei: {
        name: "BpuEI",
        seq: "CTTGAG",
        cuts: [16],
        dam: false,
        dcm: false
    },
    bsaai: {
        name: "BsaAI",
        seq: "YACGTR",
        cuts: [3],
        dam: false,
        dcm: false
    },
    bsabi: {
        name: "BsaBI",
        seq: "GATNNNNATC",
        cuts: [5],
        dam: true,
        dcm: false
    },
    bsahi: {
        name: "BsaHI",
        seq: "GRCGYC",
        cuts: [2],
        dam: false,
        dcm: true
    },
    bsai: {
        name: "BsaI",
        seq: "GGTCTC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    bsaji: {
        name: "BsaJI",
        seq: "CCNNGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsawi: {
        name: "BsaWI",
        seq: "WCCGGW",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsaxi: {
        name: "BsaXI",
        seq: "",
        cuts: [9],
        dam: false,
        dcm: false
    },
    bseri: {
        name: "BseRI",
        seq: "GAGGAG",
        cuts: [10],
        dam: false,
        dcm: false
    },
    bseyi: {
        name: "BseYI",
        seq: "CCCAGC",
        cuts: [-5],
        dam: false,
        dcm: false
    },
    bsgi: {
        name: "BsgI",
        seq: "GTGCAG",
        cuts: [16],
        dam: false,
        dcm: false
    },
    bsiei: {
        name: "BsiEI",
        seq: "CGRYCG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    bsihkai: {
        name: "BsiHKAI",
        seq: "GWGCWC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bsiwi: {
        name: "BsiWI",
        seq: "CGTACG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsli: {
        name: "BslI",
        seq: "CCNNNNNNNGG",
        cuts: [7],
        dam: false,
        dcm: true
    },
    bsmai: {
        name: "BsmAI",
        seq: "GTCTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsmbi: {
        name: "BsmBI",
        seq: "CGTCTC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    bsmfi: {
        name: "BsmFI",
        seq: "GGGAC",
        cuts: [10],
        dam: false,
        dcm: true
    },
    bsmi: {
        name: "BsmI",
        seq: "GAATGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsobi: {
        name: "BsoBI",
        seq: "CYCGRG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bspcni: {
        name: "BspCNI",
        seq: "CTCAG",
        cuts: [9],
        dam: false,
        dcm: false
    },
    bspdi: {
        name: "BspDI",
        seq: "ATCGAT",
        cuts: [2],
        dam: true,
        dcm: false
    },
    bspei: {
        name: "BspEI",
        seq: "TCCGGA",
        cuts: [1],
        dam: true,
        dcm: false
    },
    bsphi: {
        name: "BspHI",
        seq: "TCATGA",
        cuts: [1],
        dam: true,
        dcm: false
    },
    bsp1286i: {
        name: "Bsp1286I",
        seq: "GDGCHC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    bspmi: {
        name: "BspMI",
        seq: "ACCTGC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    bspqi: {
        name: "BspQI",
        seq: "GCTCTTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsrbi: {
        name: "BsrBI",
        seq: "CCGCTC",
        cuts: [-3],
        dam: false,
        dcm: false
    },
    bsrdi: {
        name: "BsrDI",
        seq: "GCAATG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsrfi: {
        name: "BsrFI",
        seq: "RCCGGY",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsrgi: {
        name: "BsrGI",
        seq: "TGTACA",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsri: {
        name: "BsrI",
        seq: "ACTGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsshii: {
        name: "BssHII",
        seq: "GCGCGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bsssi: {
        name: "BssSI",
        seq: "CACGAG",
        cuts: [-5],
        dam: false,
        dcm: false
    },
    bstbi: {
        name: "BstBI",
        seq: "TTCGAA",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bsteii: {
        name: "BstEII",
        seq: "GGTNACC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bstni: {
        name: "BstNI",
        seq: "CCWGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bstui: {
        name: "BstUI",
        seq: "CGCG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    bstyi: {
        name: "BstYI",
        seq: "RGATCY",
        cuts: [1],
        dam: false,
        dcm: false
    },
    bstz17i: {
        name: "BstZ17I",
        seq: "GTATAC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    bsu36i: {
        name: "Bsu36I",
        seq: "CCTNAGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    btgi: {
        name: "BtgI",
        seq: "CCRYGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    btgzi: {
        name: "BtgZI",
        seq: "GCGATG",
        cuts: [10],
        dam: false,
        dcm: false
    },
    btsci: {
        name: "BtsCI",
        seq: "GGATG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    btsimuti: {
        name: "BtsIMutI",
        seq: "CAGTG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    btsi: {
        name: "BtsI",
        seq: "GCAGTG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    cac8i: {
        name: "Cac8I",
        seq: "GCNNGC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    clai: {
        name: "ClaI",
        seq: "ATCGAT",
        cuts: [2],
        dam: true,
        dcm: false
    },
    cspci: {
        name: "CspCI",
        seq: "",
        cuts: [11],
        dam: false,
        dcm: false
    },
    'cviki-1': {
        name: "CviKI-1",
        seq: "RGCY",
        cuts: [2],
        dam: false,
        dcm: false
    },
    cviqi: {
        name: "CviQI",
        seq: "GTAC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    ddei: {
        name: "DdeI",
        seq: "CTNAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    dpnii: {
        name: "DpnII",
        seq: "GATC",
        cuts: [0],
        dam: true,
        dcm: false
    },
    drdi: {
        name: "DrdI",
        seq: "GACNNNNNNGTC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    eaei: {
        name: "EaeI",
        seq: "YGGCCR",
        cuts: [1],
        dam: false,
        dcm: true
    },
    eagi: {
        name: "EagI",
        seq: "CGGCCG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    eari: {
        name: "EarI",
        seq: "CTCTTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    ecii: {
        name: "EciI",
        seq: "GGCGGA",
        cuts: [11],
        dam: false,
        dcm: false
    },
    eco53ki: {
        name: "Eco53kI",
        seq: "GAGCTC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    econi: {
        name: "EcoNI",
        seq: "CCTNNNNNAGG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    ecop15i: {
        name: "EcoP15I",
        seq: "CAGCAG",
        cuts: [25],
        dam: false,
        dcm: false
    },
    fnu4hi: {
        name: "Fnu4HI",
        seq: "GCNGC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    fsei: {
        name: "FseI",
        seq: "GGCCGGCC",
        cuts: [6],
        dam: false,
        dcm: true
    },
    fspi: {
        name: "FspI",
        seq: "TGCGCA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    haeii: {
        name: "HaeII",
        seq: "RGCGCY",
        cuts: [5],
        dam: false,
        dcm: false
    },
    hinp1i: {
        name: "HinP1I",
        seq: "GCGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    hpyav: {
        name: "HpyAV",
        seq: "CCTTC",
        cuts: [6],
        dam: false,
        dcm: false
    },
    hpych4iii: {
        name: "HpyCH4III",
        seq: "ACNGT",
        cuts: [3],
        dam: false,
        dcm: false
    },
    hpych4iv: {
        name: "HpyCH4IV",
        seq: "ACGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    hpych4v: {
        name: "HpyCH4V",
        seq: "TGCA",
        cuts: [2],
        dam: false,
        dcm: false
    },
    hpy99i: {
        name: "Hpy99I",
        seq: "CGWCG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    hpy188i: {
        name: "Hpy188I",
        seq: "TCNGA",
        cuts: [3],
        dam: true,
        dcm: false
    },
    hpy166ii: {
        name: "Hpy166II",
        seq: "GTNNAC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    hpy188iii: {
        name: "Hpy188III",
        seq: "TCNNGA",
        cuts: [2],
        dam: true,
        dcm: false
    },
    kasi: {
        name: "KasI",
        seq: "GGCGCC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    mluci: {
        name: "MluCI",
        seq: "AATT",
        cuts: [0],
        dam: false,
        dcm: false
    },
    mlyi: {
        name: "MlyI",
        seq: "GAGTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    mmei: {
        name: "MmeI",
        seq: "TCCRAC",
        cuts: [20],
        dam: false,
        dcm: false
    },
    msci: {
        name: "MscI",
        seq: "TGGCCA",
        cuts: [3],
        dam: false,
        dcm: true
    },
    msei: {
        name: "MseI",
        seq: "TTAA",
        cuts: [1],
        dam: false,
        dcm: false
    },
    msli: {
        name: "MslI",
        seq: "CAYNNNNRTG",
        cuts: [5],
        dam: false,
        dcm: false
    },
    mspji: {
        name: "MspJI",
        seq: "CNNR",
        cuts: [9],
        dam: false,
        dcm: false
    },
    mwoi: {
        name: "MwoI",
        seq: "GCNNNNNNNGC",
        cuts: [7],
        dam: false,
        dcm: false
    },
    naei: {
        name: "NaeI",
        seq: "GCCGGC",
        cuts: [3],
        dam: false,
        dcm: false
    },
    nari: {
        name: "NarI",
        seq: "GGCGCC",
        cuts: [2],
        dam: false,
        dcm: false
    },
    ncii: {
        name: "NciI",
        seq: "CCSGG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    ngomiv: {
        name: "NgoMIV",
        seq: "GCCGGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    nlaiii: {
        name: "NlaIII",
        seq: "CATG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    nlaiv: {
        name: "NlaIV",
        seq: "GGNNCC",
        cuts: [3],
        dam: false,
        dcm: true
    },
    nmeaiii: {
        name: "NmeAIII",
        seq: "GCCGAG",
        cuts: [21],
        dam: false,
        dcm: false
    },
    nsii: {
        name: "NsiI",
        seq: "ATGCAT",
        cuts: [5],
        dam: false,
        dcm: false
    },
    nspi: {
        name: "NspI",
        seq: "RCATGY",
        cuts: [5],
        dam: false,
        dcm: false
    },
    paer7i: {
        name: "PaeR7I",
        seq: "CTCGAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    paqci: {
        name: "PaqCI",
        seq: "CACCTGC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    pflfi: {
        name: "PflFI",
        seq: "GACNNNGTC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    pflmi: {
        name: "PflMI",
        seq: "CCANNNNNTGG",
        cuts: [7],
        dam: false,
        dcm: true
    },
    plei: {
        name: "PleI",
        seq: "GAGTC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    pluti: {
        name: "PluTI",
        seq: "GGCGCC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    pmei: {
        name: "PmeI",
        seq: "GTTTAAAC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    pmli: {
        name: "PmlI",
        seq: "CACGTG",
        cuts: [3],
        dam: false,
        dcm: false
    },
    ppumi: {
        name: "PpuMI",
        seq: "RGGWCCY",
        cuts: [2],
        dam: false,
        dcm: true
    },
    pshai: {
        name: "PshAI",
        seq: "GACNNNNGTC",
        cuts: [5],
        dam: false,
        dcm: false
    },
    pspgi: {
        name: "PspGI",
        seq: "CCWGG",
        cuts: [0],
        dam: false,
        dcm: true
    },
    rsrii: {
        name: "RsrII",
        seq: "CGGWCCG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    sacii: {
        name: "SacII",
        seq: "CCGCGG",
        cuts: [4],
        dam: false,
        dcm: false
    },
    sapi: {
        name: "SapI",
        seq: "GCTCTTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    sau3ai: {
        name: "Sau3AI",
        seq: "GATC",
        cuts: [0],
        dam: false,
        dcm: false
    },
    sau96i: {
        name: "Sau96I",
        seq: "GGNCC",
        cuts: [1],
        dam: false,
        dcm: true
    },
    scrfi: {
        name: "ScrFI",
        seq: "CCNGG",
        cuts: [2],
        dam: false,
        dcm: true
    },
    sexai: {
        name: "SexAI",
        seq: "ACCWGGT",
        cuts: [1],
        dam: false,
        dcm: true
    },
    sfci: {
        name: "SfcI",
        seq: "CTRYAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    sfoi: {
        name: "SfoI",
        seq: "GGCGCC",
        cuts: [3],
        dam: false,
        dcm: true
    },
    sgrai: {
        name: "SgrAI",
        seq: "CRCCGGYG",
        cuts: [2],
        dam: false,
        dcm: false
    },
    smli: {
        name: "SmlI",
        seq: "CTYRAG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    snabi: {
        name: "SnaBI",
        seq: "TACGTA",
        cuts: [3],
        dam: false,
        dcm: false
    },
    spei: {
        name: "SpeI",
        seq: "ACTAGT",
        cuts: [1],
        dam: false,
        dcm: false
    },
    srfi: {
        name: "SrfI",
        seq: "GCCCGGGC",
        cuts: [4],
        dam: false,
        dcm: false
    },
    stui: {
        name: "StuI",
        seq: "AGGCCT",
        cuts: [3],
        dam: false,
        dcm: true
    },
    styd4i: {
        name: "StyD4I",
        seq: "CCNGG",
        cuts: [0],
        dam: false,
        dcm: true
    },
    styi: {
        name: "StyI",
        seq: "CCWWGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    swai: {
        name: "SwaI",
        seq: "ATTTAAAT",
        cuts: [4],
        dam: false,
        dcm: false
    },
    tfii: {
        name: "TfiI",
        seq: "GAWTC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    tsei: {
        name: "TseI",
        seq: "GCWGC",
        cuts: [1],
        dam: false,
        dcm: false
    },
    tsp45i: {
        name: "Tsp45I",
        seq: "GTSAC",
        cuts: [0],
        dam: false,
        dcm: false
    },
    tspmi: {
        name: "TspMI",
        seq: "CCCGGG",
        cuts: [1],
        dam: false,
        dcm: false
    },
    tspri: {
        name: "TspRI",
        seq: "NNCASTGNN",
        cuts: [9],
        dam: false,
        dcm: false
    },
    xcmi: {
        name: "XcmI",
        seq: "CCANNNNNNNNNTGG",
        cuts: [8],
        dam: false,
        dcm: false
    },
    xmni: {
        name: "XmnI",
        seq: "GAANNNNTTC",
        cuts: [5],
        dam: false,
        dcm: false
    }
}

REs = Object.fromEntries(
    Object.entries(REs_raw)
    .sort((a, b) => a[0].localeCompare(b[0]))
    .map(e => [e[0], {...e[1], seq_re: iupac_re(e[1].seq)}])
)

const REs_bottom = Object.fromEntries(
    Object.entries(REs)
        .filter(re => re[1].seq_re != rev_comp(re[1].seq_re))
        .map(re => [re[0], { ...re[1], seq: rev_comp(re[1].seq_re), seq_re: iupac_re(rev_comp(re[1].seq_re)) }])
)