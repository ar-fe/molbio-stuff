/**
 * Gel Electrophoresis Simulator
 * DNA restriction digest visualization
 */

// ============================================================================
// Constants and Configuration
// ============================================================================

// ============================================================================
// Constants and Configuration
// ============================================================================

const GEL_PCT_OPTIONS = Array.from({ length: 11 }, (_, i) => (i + 8) / 1000);

/**
 * Molecular weight standards library
 * Each standard contains: id, name, manufacturer, lengths, weights, and rendering coefficients
 */
const MW_STANDARDS = {
    NEB1kbPlus: {
        name: '1 kb Plus DNA Ladder',
        manufacturer: 'NEB',
        lengths: [10002, 8001, 6001, 5001, 4001, 3001, 2017, 1517, 1200, 1000, 900, 800, 700, 600, 517, 500, 400, 300, 200, 100],
        weightsRaw: [40, 40, 48, 40, 32, 120, 40, 57, 45, 122, 34, 31, 27, 23, 62, 62, 49, 37, 32, 61],
        thicknessCoefficient: 1.1,
        opacityCoefficient: 0.8,
    },
    NEB1kb: {
        name: '1 kb DNA Ladder',
        manufacturer: 'NEB',
        lengths: [10002, 8001, 6001, 5001, 4001, 3001, 2000, 1500, 1000, 517, 500],
        weightsRaw: [42, 42, 50, 42, 33, 125, 48, 36, 42, 21, 21],
        thicknessCoefficient: 1.1,
        opacityCoefficient: 0.8,
    },
    NEB100: {
        name: '100 bp DNA Ladder',
        manufacturer: 'NEB',
        lengths: [1517, 1200, 1000, 900, 800, 700, 600, 517, 500, 400, 300, 200, 100],
        weightsRaw: [45, 35, 95, 27, 24, 21, 18, 48.5, 48.5, 38, 29, 25, 48],
        thicknessCoefficient: 1.0,
        opacityCoefficient: 0.8,
    },
    NEB50: {
        name: '50 bp DNA Ladder',
        manufacturer: 'NEB',
        lengths: [1350, 916, 766, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50],
        weightsRaw: [103, 70, 58, 54, 50, 46, 42, 76, 34, 31, 27, 46, 57, 107, 46, 69, 84],
        thicknessCoefficient: 0.9,
        opacityCoefficient: 0.8,
    },
    generulerMix: {
        name: 'GeneRuler DNA Ladder Mix',
        manufacturer: 'Thermo',
        lengths: [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        weightsRaw: [18, 18, 18, 18, 18, 18, 60, 16, 16, 16, 16, 60, 17, 17, 17, 17, 60, 20, 20, 20, 20],
        thicknessCoefficient: 1.0,
        opacityCoefficient: 0.8,
    },
    generuler1kbPlus: {
        name: 'GeneRuler 1 kb Plus',
        manufacturer: 'Thermo',
        lengths: [20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300, 200, 75],
        weightsRaw: [20, 20, 20, 75, 20, 20, 20, 80, 25, 25, 75, 25, 25, 25, 25],
        thicknessCoefficient: 1.0,
        opacityCoefficient: 0.8,
    },
    generuler1kb: {
        name: 'GeneRuler 1 kb',
        manufacturer: 'Thermo',
        lengths: [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 750, 500, 250],
        weightsRaw: [30, 30, 70, 30, 30, 30, 70, 25, 25, 25, 60, 25, 25, 25],
        thicknessCoefficient: 1.0,
        opacityCoefficient: 0.8,
    },
    generuler100plus: {
        name: 'GeneRuler 100 bp Plus',
        manufacturer: 'Thermo',
        lengths: [3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        weightsRaw: [28, 28, 28, 28, 80, 27, 27, 27, 27, 80, 30, 30, 30, 30],
        thicknessCoefficient: 1.0,
        opacityCoefficient: 0.8,
    },
    generuler100: {
        name: 'GeneRuler 100 bp',
        manufacturer: 'Thermo',
        lengths: [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        weightsRaw: [45, 45, 45, 45, 45, 115, 40, 40, 40, 40],
        thicknessCoefficient: 1.0,
        opacityCoefficient: 0.8,
    },
    generuler50: {
        name: 'GeneRuler 50 bp',
        manufacturer: 'Thermo',
        lengths: [1000, 900, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50],
        weightsRaw: [30, 30, 30, 30, 30, 75, 30, 30, 75, 35, 35, 35, 35],
        thicknessCoefficient: 1.0,
        opacityCoefficient: 0.8,
    },
    evrogen1k: {
        name: 'DNA Ladder 1 kb',
        manufacturer: 'Evrogen',
        lengths: [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250],
        weightsRaw: [10, 10, 10, 10, 10, 20, 10, 10, 10, 20, 10, 10, 10],
        thicknessCoefficient: 1.0,
        opacityCoefficient: 0.8,
    },
    evrogen100: {
        name: 'DNA Ladder 100 bp',
        manufacturer: 'Evrogen',
        lengths: [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        weightsRaw: [7, 5, 5, 20, 5, 15, 10, 15, 9, 9],
        thicknessCoefficient: 1.2,
        opacityCoefficient: 0.7,
    },
    evrogen50: {
        name: 'DNA Ladder 50+ bp',
        manufacturer: 'Evrogen',
        lengths: [700, 500, 400, 350, 300, 250, 200, 150, 100, 50],
        weightsRaw: [30, 30, 15, 15, 30, 15, 20, 15, 15, 15],
        thicknessCoefficient: 0.9,
        opacityCoefficient: 0.8,
    },
    blmSkyHigh: {
        name: 'Sky-High',
        manufacturer: 'Biolabmix',
        lengths: [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250],
        weightsRaw: [75, 75, 75, 75, 100, 200, 70, 60, 50, 100, 30, 40, 50],
        thicknessCoefficient: 1.25,
        opacityCoefficient: 1,
    },
    blmStep100long: {
        name: 'Step100 Long',
        manufacturer: 'Biolabmix',
        lengths: [3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        weightsRaw: [179, 118, 179, 70, 71, 64, 57, 50, 43, 71, 29, 20, 29, 20],
        thicknessCoefficient: 1.25,
        opacityCoefficient: 0.8,
    },
    blmStep100: {
        name: 'Step100',
        manufacturer: 'Biolabmix',
        lengths: [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
        weightsRaw: [156, 140, 126, 110, 94, 156, 62, 46, 62, 46],
        thicknessCoefficient: 0.8,
        opacityCoefficient: 0.9,
    },
    blmStep100_50: {
        name: 'Step100+50',
        manufacturer: 'Biolabmix',
        lengths: [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 50],
        weightsRaw: [150, 130, 120, 100, 90, 150, 60, 50, 60, 50, 40],
        thicknessCoefficient: 0.8,
        opacityCoefficient: 0.9,
    },
    blmStep50plus: {
        name: 'Step50 plus',
        manufacturer: 'Biolabmix',
        lengths: [1500, 1000, 700, 600, 500, 400, 350, 300, 250, 200, 150, 100, 50],
        weightsRaw: [200, 135, 95, 80, 135, 50, 45, 40, 35, 80, 40, 40, 25],
        thicknessCoefficient: 1.25,
        opacityCoefficient: 0.9,
    },
};

/**
 * Get MW standard with computed weights (normalized by length)
 */
const getMWStandard = (standardId) => {
    const standard = MW_STANDARDS[standardId];
    if (!standard) return null;

    return {
        ...standard,
        weights: standard.weightsRaw.map((w, i) => w / standard.lengths[i]),
    };
};

const DEFAULT_MW_STANDARD = 'evrogen1k';
let currentMWStandard = DEFAULT_MW_STANDARD;

// ============================================================================
// State
// ============================================================================

let digestFragmentLengths = [];
let digestFragmentIntensities = [];
let circularProbability = 0;

// ============================================================================
// Intensity Calibration Parameters
// ============================================================================

/**
 * Intensity rendering model configuration
 * 
 * Physics basis:
 * - Larger fragments stain better (log correction)
 * - Band brightness depends on DNA quantity and stain accessibility
 * - Different lane types (marker vs sample) often loaded differently
 */
const INTENSITY_MODEL = {
    // Size-dependent correction: how fragment length affects staining intensity
    sizeFactorMode: 'none',  // 'none' | 'sqrt' | 'log'
    
    // Brightness scaling per lane type
    lanes: {
        mwStandard: {
            maxBandThickness: 0.25,    // maximum band thickness in %
            minBandThickness: 0.08,    // minimum band thickness in %
            maxOpacity: 0.9,           // opacity at max intensity
            minOpacity: 0.35,          // opacity at min intensity
        },
        digestedDna: {
            maxBandThickness: 0.4,     // maximum band thickness in %
            minBandThickness: 0.12,    // minimum band thickness in %
            maxOpacity: 0.75,          // opacity at max intensity
            minOpacity: 0.05,          // opacity at min intensity
        },
    },
    
    // SVG rendering parameters
    intensityNormalizationFactor: 2.5,          // scaling factor for migration distance
    maxTotalOpacity: 0.95,                      // clip individual opacity to this value
    maxCumulativeOpacity: 5,                    // clip sum of all opacities to this value
};

/**
 * Apply size-dependent correction to fragment intensity
 * Larger fragments tend to stain better, smaller fragments are harder to see
 */
const applySizeCorrectionFactor = (fragmentLength, rawIntensity) => {
    switch (INTENSITY_MODEL.sizeFactorMode) {
        case 'sqrt':
            return rawIntensity * Math.sqrt(fragmentLength / 1000);  // normalize to 1kb
        case 'log':
            return rawIntensity * (0.5 + 0.5 * Math.log10(fragmentLength) / 4);  // log-corrected, clipped to [0.5, 1.5]
        case 'none':
        default:
            return rawIntensity;
    }
};

// ============================================================================
// Migration Distance Calculations
// ============================================================================

/**
 * Calculate migration distance for a fragment of given length in gel at given percentage
 */
const calculateMigrationDistance = (fragmentLength, gelPercentage) => {
    const L = 0.375 + gelPercentage;
    const x0 = 3.2 - 17 * gelPercentage;
    const k = -1.9 - 42 * gelPercentage;
    const b = 0.25 - 8.8 * gelPercentage;
    const x = Math.log10(fragmentLength);
    return L / (1 + Math.exp(-k * (x - x0))) + b;
};

const calculateInfinityMigrationDistance = (gelPercentage) => 0.25 - 8.8 * gelPercentage;
const calculateZeroMigrationDistance = (gelPercentage) => 0.625 - 7.8 * gelPercentage;

/**
 * Calculate normalized migration distance (0 to 1 scale)
 */
const calculateNormalizedMigrationDistance = (fragmentLength, gelPercentage) => {
    const infinityDistance = calculateInfinityMigrationDistance(gelPercentage);
    const zeroDistance = calculateZeroMigrationDistance(gelPercentage);
    return (calculateMigrationDistance(fragmentLength, gelPercentage) - infinityDistance) /
           (zeroDistance - infinityDistance);
};

// ============================================================================
// SVG Band Generation
// ============================================================================

/**
 * Generate SVG rectangles for gel bands
 * @param {number} gelPercentage - Gel concentration percentage
 * @param {number} migrationTime - Time parameter for migration
 * @param {number[]} fragmentLengths - DNA fragment lengths
 * @param {number[]} fragmentIntensities - Band intensities
 * @param {number} intensityCoefficient - Scaling coefficient for intensities
 * @param {boolean} isMwStandard - true if rendering MW marker
 * @param {object} mwStandardConfig - MW standard config with rendering coefficients
 */
const generateBandsSvg = (
    gelPercentage,
    migrationTime,
    fragmentLengths,
    fragmentIntensities = null,
    intensityCoefficient = 1,
    isMwStandard = false,
    mwStandardConfig = null,
) => {
    const laneConfig = isMwStandard ? INTENSITY_MODEL.lanes.mwStandard : INTENSITY_MODEL.lanes.digestedDna;
    const weights = fragmentIntensities || fragmentLengths.map(() => 1);
    const intensities = weights.map((w, i) => w * fragmentLengths[i]);
    
    // Apply size correction
    const sizeAdjustedIntensities = intensities.map((intensity, i) =>
        intensity * applySizeCorrectionFactor(fragmentLengths[i], 1),
    );
    
    // Calculate band thicknesses: proportional to intensity
    // For MW: use absolute intensities to show 2x bands as 2x thick
    // For sample: use relative intensities
    const maxIntensity = Math.max(...sizeAdjustedIntensities);
    const bandThicknesses = sizeAdjustedIntensities.map((i, idx) => {
        let ratio;
        if (isMwStandard) {
            // Absolute: preserve intensity differences (2x intensity = 2x thick)
            ratio = i / maxIntensity;
        } else {
            // Relative: normalize to show proportion
            const totalIntensity = sizeAdjustedIntensities.reduce((a, b) => a + b, 0);
            ratio = (i / totalIntensity) * intensityCoefficient;
            ratio = Math.min(ratio * 2, 1);  // scale up for visibility, but cap at 1
        }
        let thickness = laneConfig.minBandThickness + ratio * (laneConfig.maxBandThickness - laneConfig.minBandThickness);
        // Apply MW standard-specific thickness coefficient
        if (isMwStandard && mwStandardConfig) {
            thickness *= mwStandardConfig.thicknessCoefficient;
        }
        return thickness;
    });

    // Calculate band opacities: proportional to intensity
    let bandOpacities = sizeAdjustedIntensities.map(i => {
        const ratio = i / maxIntensity;
        let opacity = laneConfig.minOpacity + ratio * (laneConfig.maxOpacity - laneConfig.minOpacity);
        // Apply MW standard-specific opacity coefficient
        if (isMwStandard && mwStandardConfig) {
            opacity *= mwStandardConfig.opacityCoefficient;
        }
        // Clip individual opacity to reasonable limit
        return Math.min(opacity, INTENSITY_MODEL.maxTotalOpacity);
    });
    
    // Normalize cumulative opacity if sum exceeds limit
    const totalOpacity = bandOpacities.reduce((a, b) => a + b, 0);
    if (totalOpacity > INTENSITY_MODEL.maxCumulativeOpacity && !isMwStandard) {
        const scaleFactor = INTENSITY_MODEL.maxCumulativeOpacity / totalOpacity;
        bandOpacities = bandOpacities.map(o => o * scaleFactor);
    }
    
    const migrationDistances = fragmentLengths.map(length =>
        calculateNormalizedMigrationDistance(length, gelPercentage),
    );

    return fragmentLengths
        .map(
            (length, index) =>
                `<rect width="2" height="${bandThicknesses[index]}" y="-${
                    migrationDistances[index] * INTENSITY_MODEL.intensityNormalizationFactor * migrationTime
                }%" fill="rgba(0, 0, 0, ${bandOpacities[index]})" />`,
        )
        .join('\n');
};

// ============================================================================
// Gel Redraw
// ============================================================================

const redrawGel = () => {
    const gelPercentage = +document.getElementById('gel-pct').value;
    const migrationTime = document.getElementById('gel-time').value;
    const mwStandard = getMWStandard(currentMWStandard);

    const mwMarkerLaneElement = document.getElementById('left-lane');
    mwMarkerLaneElement.innerHTML = generateBandsSvg(
        gelPercentage,
        migrationTime,
        mwStandard.lengths,
        mwStandard.weights,
        3,
        true,  // isMwStandard flag
        mwStandard,  // pass MW standard config with coefficients
    );

    const digestLaneElement = document.getElementById('right-lane');
    digestLaneElement.innerHTML = generateBandsSvg(
        gelPercentage,
        migrationTime,
        digestFragmentLengths,
        digestFragmentIntensities,
        (1 - circularProbability) / 2.5,
        false,  // isMwStandard flag
        null,   // no config for sample
    );

    document.getElementById('question-mark').style.fill = `rgba(0, 0, 0, ${circularProbability})`;
};

const digestAndRedraw = () => {
    digestDNA();
    redrawGel();
};

// ============================================================================
// Restriction Enzyme List Management
// ============================================================================

/**
 * Get restriction enzymes with their efficiencies from UI
 */
const getSelectedEnzymesWithEfficiency = () =>
    Array.from(document.querySelectorAll('.re-box input.valid-re')).map(el => ({
        enzymeName: el.value.toLowerCase(),
        efficiency: +el.parentElement.parentElement.querySelector('.re-eff').value,
    }));

const getSelectedEnzymeNames = () =>
    Array.from(new Set(getSelectedEnzymesWithEfficiency().map(e => e.enzymeName)));

/**
 * Update visibility of enzyme dropdown options based on user input
 */
const updateEnzymeDropdownOptions = () => {
    const availableEnzymes = Object.keys(REs);
    const selectedEnzymes = getSelectedEnzymeNames();

    document.querySelectorAll('.re-box').forEach(box => {
        const inputElement = box.querySelector('input');
        const dropdownElement = box.querySelector('.re-list');
        const userInput = inputElement.value.toLowerCase();
        const matchingEnzymes = availableEnzymes.filter(enzyme => enzyme.includes(userInput));

        dropdownElement.childNodes.forEach(option => {
            option.hidden = !(
                matchingEnzymes.includes(option.value) &&
                (!selectedEnzymes.includes(option.value) || option.value === userInput)
            );
        });
    });
};

// ============================================================================
// Sequence Utilities
// ============================================================================

const getSequenceFromInput = () =>
    document.getElementById('seq').value.replace(/[^ACGTacgt]/g, '').toUpperCase();

const clearSequenceInput = () => {
    const sequenceElement = document.getElementById('seq');
    sequenceElement.value = '';
    sequenceElement.dispatchEvent(new Event('input'));
};

const getReverseComplement = (dnaSequence) => {
    const complementMap = { A: 'T', C: 'G', G: 'C', T: 'A' };
    return dnaSequence
        .split('')
        .map(base => complementMap[base])
        .reverse()
        .join('');
};

// ============================================================================
// Combinatorial Utilities
// ============================================================================

/**
 * Generate all subsets of an array (power set)
 */
const generatePowerSet = (array) => {
    const subsets = [[]];
    for (const element of array) {
        const currentLength = subsets.length;
        for (let i = 0; i < currentLength; i++) {
            subsets.push([...subsets[i], element]);
        }
    }
    return subsets;
};

/**
 * Combine two sets with probabilities
 */
const combineWithProbabilities = (set1, set2) => {
    if (!set1.length) return set2;
    if (!set2.length) return set1;

    const result = [];
    set1.forEach(item1 => {
        set2.forEach(item2 => {
            const mergedEnds = [...item1.ends, ...item2.ends].sort((a, b) => a - b);
            result.push({ ends: mergedEnds, probability: item1.probability * item2.probability });
        });
    });
    return result;
};

/**
 * Compute cartesian product with probabilities for multiple arrays
 */
const computeCartesianProductWithProbabilities = (arrays) =>
    arrays.length === 0 ? [] : arrays.reduce((acc, curr) => combineWithProbabilities(acc, curr), []);

// ============================================================================
// DNA Cutting
// ============================================================================

/**
 * Cut DNA at specified positions for all selected restriction enzymes
 */
const cutDNAByEnzymes = (sequence, selectedEnzymes, modifiedDamSequence, modifiedDcmSequence) => {
    const restrictionCuts = { LINEAR: [{ ends: [0, sequence.length], probability: 1 }] };

    for (const enzyme of selectedEnzymes) {
        const enzymeConfig = REs[enzyme.enzymeName];
        const efficiency = enzyme.efficiency;

        let workingSequence = sequence;
        if (enzymeConfig.dam && document.getElementById('dam').checked) workingSequence = modifiedDamSequence;
        if (enzymeConfig.dcm && document.getElementById('dcm').checked) workingSequence = modifiedDcmSequence;

        const topStrandCutPositions = [...workingSequence.matchAll(enzymeConfig.seq_re, 'g')].flatMap(match =>
            enzymeConfig.cuts.map(cutOffset => match.index + cutOffset),
        );

        let bottomStrandCutPositions = [];
        if (enzyme.enzymeName in REs_bottom) {
            const enzymeConfigBottom = REs_bottom[enzyme.enzymeName];
            bottomStrandCutPositions = [...workingSequence.matchAll(enzymeConfigBottom.seq_re, 'g')].flatMap(match =>
                enzymeConfigBottom.cuts.map(cutOffset => match.index - cutOffset),
            );
        }

        const allCutPositions = Array.from(new Set([...topStrandCutPositions, ...bottomStrandCutPositions])).sort(
            (a, b) => a - b,
        );
        allCutPositions.unshift(0);
        allCutPositions.push(sequence.length);

        const cutSubsets = efficiency === 1 ? [allCutPositions] : generatePowerSet(allCutPositions);
        restrictionCuts[enzyme.enzymeName] = cutSubsets.map(cuts => ({
            ends: cuts,
            probability:
                Math.pow(efficiency, cuts.length) * Math.pow(1 - efficiency, allCutPositions.length - cuts.length),
        }));
    }

    const allPossibleFragments = computeCartesianProductWithProbabilities(Object.values(restrictionCuts));
    if (allPossibleFragments.length === 0) return allPossibleFragments;

    return allPossibleFragments.map(fragmentSet => {
        const fragmentPairs = fragmentSet.ends.slice(1).map((end, i) => [fragmentSet.ends[i], end]);
        return {
            lengths: fragmentPairs.map(pair => pair[1] - pair[0]).filter(l => l > 0),
            probability: fragmentSet.probability,
        };
    });
};

/**
 * Main digestion function
 */
const digestDNA = () => {
    const sequence = getSequenceFromInput();
    const sequenceLength = sequence.length;

    if (sequenceLength === 0) {
        digestFragmentLengths = [];
        digestFragmentIntensities = [];
        circularProbability = 0;
        return;
    }

    const selectedEnzymes = getSelectedEnzymesWithEfficiency();
    const isCircular = document.getElementById('circ').checked;
    circularProbability = 0;

    const connectCircularFragments = (lengths, probability) => {
        if (lengths.length === 0) return [];
        if (lengths.length === 1) circularProbability += probability;

        const connectedLength = lengths.at(0) + lengths.at(-1);
        const middleFragments = lengths.slice(1, -1);
        middleFragments.push(connectedLength);
        return middleFragments;
    };

    const damModifiedSequence = sequence.repeat(2).replace(/GATC/g, 'GmmC').slice(5, sequenceLength + 5);
    const dcmModifiedSequence = sequence.repeat(2).replace(/CC([AT])GG/g, 'Cm$1mG').slice(5, sequenceLength + 5);

    const fragmentPieces = cutDNAByEnzymes(sequence, selectedEnzymes, damModifiedSequence, dcmModifiedSequence).map(
        piece => ({
            lengths: isCircular ? connectCircularFragments(piece.lengths, piece.probability) : piece.lengths,
            probability: piece.probability,
        }),
    );

    const fragmentLengthProbabilityMap = new Map();
    fragmentPieces.forEach(piece => {
        piece.lengths.forEach(length => {
            if (length > 0 && length <= sequenceLength) {
                fragmentLengthProbabilityMap.set(
                    length,
                    (fragmentLengthProbabilityMap.get(length) || 0) + piece.probability,
                );
            }
        });
    });

    digestFragmentLengths = Array.from(fragmentLengthProbabilityMap.keys());
    digestFragmentIntensities = Array.from(fragmentLengthProbabilityMap.values());
    const totalIntensity = digestFragmentIntensities.reduce((a, b) => a + b, 0);

    if (digestFragmentIntensities.length === 0) {
        circularProbability = 1;
    } else {
        digestFragmentIntensities = digestFragmentIntensities.map(w => w / totalIntensity);
    }
};

// ============================================================================
// Helper Functions
// ============================================================================

const clearEnzymeInput = (clearButton) => {
    const inputElement = clearButton.parentElement.querySelector('input');
    inputElement.value = '';
    const efficiencySelect = clearButton.parentElement.parentElement.querySelector('.re-eff');
    efficiencySelect.value = 1;
    inputElement.dispatchEvent(new Event('input'));
};

// ============================================================================
// Event Handler Initialization
// ============================================================================

const initializeGelPercentageSelect = () => {
    const gelPercentageElement = document.getElementById('gel-pct');
    if (gelPercentageElement) {
        gelPercentageElement.innerHTML = GEL_PCT_OPTIONS
            .map(pct => `<option value="${pct}">${(pct * 100).toFixed(1).replace('.', ',')} %</option>`)
            .join('');
        gelPercentageElement.value = '0.01';
        gelPercentageElement.addEventListener('input', redrawGel);
    }
};

const initializeEnzymeDropdowns = () => {
    document.querySelectorAll('.re-box .re-list').forEach(dropdownElement => {
        dropdownElement.innerHTML = Object.entries(REs)
            .map(([enzymeKey, enzymeData]) => `<option value="${enzymeKey}">${enzymeData.name}</option>`)
            .join('\n');
    });
};

const initializeEfficiencySelects = () => {
    const efficiencyOptions = [1.0, 0.9, 0.8, 0.5, 0.25, 0.1];
    document.querySelectorAll('.re-box .re-eff').forEach(selectElement => {
        selectElement.innerHTML = efficiencyOptions
            .map(eff => `<option value="${eff}">${eff * 100} %</option>`)
            .join('\n');
        selectElement.addEventListener('input', digestAndRedraw);
        // selectElement.addEventListener('input', e => {
        //     const el = e.currentTarget
        //     el.style.backgroundColor = `rgba(255, 0, 0, ${(1 - el.value) / 2})`
        // });
    });
};

const attachEnzymeDropdownBlurHandlers = () => {
    document.querySelectorAll('.re-box .re-list').forEach(dropdownElement => {
        dropdownElement.addEventListener('blur', () => {
            dropdownElement.selectedIndex = -1;
        });
    });
};

const attachEnzymeDropdownClickHandlers = () => {
    document.querySelectorAll('.re-box .re-list').forEach(dropdownElement => {
        dropdownElement.addEventListener('click', e => {
            const inputElement = e.currentTarget.parentElement.querySelector('input');
            if (dropdownElement.selectedIndex >= 0) {
                inputElement.value = dropdownElement.options[dropdownElement.selectedIndex].value;
                inputElement.dispatchEvent(new Event('input'));
            }
        });
    });
};

const attachEnzymeInputHandlers = () => {
    document.querySelectorAll('.re-box input').forEach(inputElement => {
        inputElement.addEventListener('input', e => {
            const currentInput = e.currentTarget;
            const efficiencySelect = currentInput.parentElement.parentElement.querySelector('.re-eff');
            const userInput = currentInput.value.toLowerCase();

            if (userInput === '') {
                currentInput.classList.remove('valid-re', 'invalid-re');
                efficiencySelect.disabled = true;
                updateEnzymeDropdownOptions();
                digestAndRedraw();
                return;
            }

            const availableEnzymes = Object.keys(REs);
            if (availableEnzymes.includes(userInput)) {
                currentInput.classList.add('valid-re');
                currentInput.classList.remove('invalid-re');
                efficiencySelect.disabled = false;
            } else {
                currentInput.classList.add('invalid-re');
                currentInput.classList.remove('valid-re');
                efficiencySelect.disabled = true;
            }

            updateEnzymeDropdownOptions();
            digestAndRedraw();
        });

        // Arrow down autofill: fill with first visible option in dropdown
        inputElement.addEventListener('keydown', e => {
            if (e.key === 'ArrowDown') {
                e.preventDefault();
                const dropdownElement = inputElement.parentElement.parentElement.querySelector('.re-list');
                const visibleOptions = Array.from(dropdownElement.options).filter(opt => !opt.hidden);
                if (visibleOptions.length > 0) {
                    inputElement.value = visibleOptions[0].value;
                    inputElement.dispatchEvent(new Event('input'));
                }
            }
        });
    });
};

const attachCircularLinearToggleHandlers = () => {
    document.querySelectorAll('input[name="circ-lin"]').forEach(radioElement => {
        radioElement.addEventListener('input', e => {
            const selectedRadio = e.currentTarget;
            const damCheckbox = document.getElementById('dam');
            const dcmCheckbox = document.getElementById('dcm');

            if (selectedRadio.id === 'lin') {
                damCheckbox.checked = false;
                dcmCheckbox.checked = false;
            } else if (selectedRadio.id === 'circ') {
                damCheckbox.checked = true;
                dcmCheckbox.checked = true;
            }
            digestAndRedraw();
        });
    });
};

const attachSequenceInputHandler = () => {
    document.getElementById('seq').addEventListener('input', digestAndRedraw);
};

const attachModificationCheckboxHandlers = () => {
    document.getElementById('dam').addEventListener('input', digestAndRedraw);
    document.getElementById('dcm').addEventListener('input', digestAndRedraw);
};

const attachClearSequenceHandler = () => {
    document.getElementById('clear-seq-btn').addEventListener('click', clearSequenceInput);
};

const attachClearEnzymeHandlers = () => {
    document.querySelectorAll('.clear-enzyme-btn').forEach(buttonElement => {
        buttonElement.addEventListener('click', e => {
            clearEnzymeInput(e.currentTarget);
        });
    });
};

const attachGelTimeSliderHandler = () => {
    document.getElementById('gel-time').addEventListener('input', redrawGel);
};

const initializeMWStandardSelect = () => {
    const mwStandardSelect = document.getElementById('mw-standard');
    if (!mwStandardSelect) return;

    let manufacturers = []
    for (let standard in MW_STANDARDS) {
        if (!manufacturers.includes(MW_STANDARDS[standard].manufacturer)) {
            manufacturers.push(MW_STANDARDS[standard].manufacturer)
        }
    }

    let standardGroups = []
    for (let manufacturer of manufacturers) {
        let standards = Object.fromEntries(Object.entries(MW_STANDARDS).filter(mw => mw[1].manufacturer == manufacturer))
        let options = Object.entries(standards)
            .map(([id, standard]) => `<option value="${id}">${standard.name}</option>`)
            .join('');
        let group = `<optgroup label="${manufacturer}">${options}</optgroup>`
        standardGroups.push(group)
    }

    // mwStandardSelect.innerHTML = Object.entries(MW_STANDARDS)
    //     .map(([id, standard]) => `<option value="${id}">${standard.name}</option>`)
    //     .join('');
    mwStandardSelect.innerHTML = standardGroups.join('');
    
    mwStandardSelect.value = DEFAULT_MW_STANDARD;
    mwStandardSelect.addEventListener('change', e => {
        currentMWStandard = e.currentTarget.value;
        redrawGel();
    });
};

// ============================================================================
// Application Startup
// ============================================================================

document.addEventListener('DOMContentLoaded', () => {
    initializeGelPercentageSelect();
    initializeEnzymeDropdowns();
    attachEnzymeDropdownBlurHandlers();
    attachEnzymeDropdownClickHandlers();
    initializeEfficiencySelects();
    attachEnzymeInputHandlers();
    attachCircularLinearToggleHandlers();
    attachSequenceInputHandler();
    attachModificationCheckboxHandlers();
    attachClearSequenceHandler();
    attachClearEnzymeHandlers();
    attachGelTimeSliderHandler();
    initializeMWStandardSelect();
    redrawGel();
});