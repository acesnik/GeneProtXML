# module name: genemodel
# main program: geneprotxml

import re
import math

#IUPAC letters
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
amino_acids_1to3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys',
    'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val',
    'W': 'Trp', 'Y': 'Tyr'
}
amino_acids_3to1 = dict((x[1], x[0]) for x in amino_acids_1to3.items())

ambiguous_dna_letters = "GATCRYWSMKHBVDN"
unambiguous_dna_letters = "GATC"
ambiguous_dna_values = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "M": "AC", "R": "AG", "W": "AT", "S": "CG", "Y": "CT", "K": "GT",
    "V": "ACG", "H": "ACT", "D": "AGT", "B": "CGT", "X": "GATC", "N": "GATC"
}
ambiguous_dna_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N"
}
unambiguous_dna_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}
standard_start_codons = ['ATG']
extended_start_codons = ['TTG', 'CTG', 'ATG']
standard_stop_codons = ['TAA', 'TAG', 'TGA']
mitochon_start_codons = ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
mitochon_stop_codons = ['TAA', 'TAG', 'AGA', 'AGG']
standard_base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
standard_base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
standard_base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
standard_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
mitochon_acids = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"

standard_code = {standard_base1[i] + standard_base2[i] + standard_base3[i]: standard_acids[i] for i in range(len(standard_acids))}
mitochon_code = {standard_base1[i] + standard_base2[i] + standard_base3[i]: mitochon_acids[i] for i in range(len(mitochon_acids))}


def translate(dna_seq):
    n = len(dna_seq)
    translation = ''
    for i in range(0, n - n % 3, 3):
        codon = dna_seq[i: i+3]
        if codon not in standard_code: return None
        else: translation += standard_code[codon]
    return translation


def reverse_complement(dna_seq):
    rc = ''
    for letter in dna_seq[::-1]:
        if letter not in unambiguous_dna_complement: return None
        else: rc += unambiguous_dna_complement[letter]
    return rc


#Gene Model structure
class Chromosome:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.genes = []
        self.amino_acid_sequences = []

    def __len__(self): return len(self.sequence)

    def __str__(self): return self.sequence

    def __contains__(self, gene_name): return gene_name in [gene.name for gene in self.genes]


class ChromSegment:
    def __init__(self, id, chrom, strand, start, end, name, biotype):
        self.name = name
        self.id = id
        self.biotype = biotype
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end

    def __len__(self): return self.end - self.start

    def __str__(self): return self.chrom.sequence[self.start: self.end + 1]

    def is_before(self, segment): return self.start < segment.start and self.end < segment.end and self.end < segment.start

    def is_after(self, segment): return self.start > segment.start and self.end > segment.end and self.start > segment.end

    def overlaps(self, segment): return not self.is_before(segment) and not self.is_after(segment)

    def equals(self, segment): return self.chrom == segment.chrom and self.start == segment.start and self.end == segment.end

    def bed_text(self, score, blockCount, blockSizes, blockStarts):
        return '\t'.join([str(self.chrom.name), str(self.start), str(self.end), self.id,
                  str(int(score)), self.strand, str(self.start), str(self.end), '255,0,0',
                        str(blockCount), ','.join([str(blockSize) for blockSize in blockSizes]),
                        ','.join([str(blockStart) for blockStart in blockStarts])])


#Variation constructs
class SequenceVariant(ChromSegment):
    def __init__(self, chrom, position, id, reference, alternate, qual, allele_frequency, depth):
        length = math.max(len(reference), len(alternate))
        ChromSegment.__init__(self, id, chrom, '+', position, position + length - 1, None, None)
        self.qual = qual
        self.reference = reference
        self.alternate = alternate
        self.allele_frequency = float(allele_frequency)
        self.depth = depth

    def __str__(self):
        return ':'.join([self.chrom.name, str(self.start), str(self.end), self.reference, self.alternate, str(self.allele_frequency)])


class SNV(SequenceVariant):
    def __init__(self, chrom, position, id, reference, alternate, qual, allele_frequency, depth):
        SequenceVariant.__init__(self, chrom, position, id, reference, alternate, qual, allele_frequency, depth)

    
    def is_missense(self):
        pass

    # TODO: replace this with a check against the gene model instead of taking it from the snpeff annotations
    def parse_aa_change(self, aa_change):
        aa_abbrev_dict = amino_acids_3to1
        aa_change_regex = '([A-Z])(\d+)([A-Z])'  # G528R
        aa_hgvs_regex = 'p\.([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])(/c\.(\d+)([ACGTN])>([ACGTN]))'  # p.Gly528Arg/c.1582G>C
        aa_pos = None  # 1-based position
        ref_aa, alt_aa = '_', '_'
        m = re.match(aa_change_regex,
                     aa_change)  # parse aa_change, and get AA change position and alternate Animo Acid
        if m:
            aa_pos = int(m.groups()[1])
            ref_aa = m.groups()[0]
            alt_aa = m.groups()[2]
        else:
            m = re.match(aa_hgvs_regex, aa_change)
            if m:
                aa_pos = int(m.groups()[1])
                ref_aa = aa_abbrev_dict[m.groups()[0]]
                alt_aa = aa_abbrev_dict[m.groups()[2]]
        return aa_pos, ref_aa, alt_aa


class Indel(ChromSegment):
    def __init__(self, chrom, position, id, reference, alternate, qual, allele_frequency, depth):
        SequenceVariant.__init__(self, chrom, position, id, reference, alternate, qual, allele_frequency, depth)


class Insertion(Indel):
    def __init__(self, chrom, position, id, reference, alternate, qual, allele_frequency, depth):
        Indel.__init__(self, chrom, position, id, reference, alternate, qual, allele_frequency, depth)


class Deletion(Indel):
    def __init__(self, chrom, position, id, reference, alternate, qual, allele_frequency, depth):
        Indel.__init__(self, chrom, position, id, reference, alternate, qual, allele_frequency, depth)


class LocalHaplotype(ChromSegment):
    def __init__(self, seqvar_list):
        if not seqvar_list or False in [isinstance(seqvar, SequenceVariant) for seqvar in seqvar_list]: return None
        self.seqvar_list = seqvar_list
        start = sorted(self.seqvar_list, key=lambda x: x.start)[0].start
        end = sorted(self.seqvar_list, key=lambda x: x.end)[0].end
        chroms = [seqvar.chrom for seqvar in self.seqvar_list]
        if len(chroms) > 1:
            print "Error: Haplotypes across multiple chromosomes not expected in object LocalHaplotype."
            exit(2)
        ChromSegment.__init__(self, None, chroms[0], '+', start, end, None, None)


#Amino acid sequences
class Modification:
    def __init__(self, description, accession, feature_type, position, target_aas, monoisotopic_mass_shift, average_mass_shift):
        self.description = description
        self.accession = accession
        self.feature_type = feature_type
        self.position = position
        self.target_aas = target_aas
        self.monoisotopic_mass_shift = monoisotopic_mass_shift
        self.average_mass_shift = average_mass_shift

    def __str__(self):
        return "description=" + self.description + " accession=" + self.accession + " feature_type=" + \
            self.feature_type + " monoisotopic_mass=" + self.monoisotopic_mass_shift


class AminoAcidSequence:
    def __init__(self, name, sequence='', transcript=None, trans_start=None, trans_end=None):
        self.name = name
        self.original_sequence = sequence
        self.sequence = sequence
        self.modifications = {} # location : modification
        self.sequence_variants = {} # location : sequence_variant
        self.other_features = [] # such as propeptides
        self.transcript = transcript
        self.trans_start = trans_start
        self.trans_end = trans_end

    def __len__(self): return len(self.sequence)

    def __str__(self): return self.sequence

    def get_seq(self): return self.sequence

    def tryptic_digestion(self):
        tryptic_peps = []
        k_fragments = self.sequence.split('K')
        for i, k_fragment in enumerate(k_fragments):
            if i != len(k_fragments) - 1:
                kr_fragments = (k_fragment + 'K').split('R')
            else:
                kr_fragments = k_fragment.split('R')
            for j, kr_fragment in enumerate(kr_fragments):
                if j != len(kr_fragments) - 1:
                    tryptic_peps.append(kr_fragment + 'R')
                else:
                    tryptic_peps.append(kr_fragment)
        while tryptic_peps and not tryptic_peps[-1]: tryptic_peps.pop()  # remove empty elements at end of list

        tryp_start, tryp_end = self.trans_start, self.trans_start
        for i, tryptic_pep in enumerate(tryptic_peps):
            tryp_end += len(tryptic_pep) * 3
            yield AminoAcidSequence(self.name + "_" + str(i), tryptic_pep, self.transcript, tryp_start, tryp_end)
            tryp_start += len(tryptic_pep) * 3

    def trim(self, start, end):
        self.sequence = self.sequence[start:end]
        self.trans_start += start * 3
        self.trans_end -= (len(self) - end) * 3

    def methionine_cleavage(self):
        if str(self)[0] == 'M': self.trim(1, len(self))
        self.trans_start += 3

    def bed_text(self):
        pass

    def lxml_element(self):
        pass





class Gene(ChromSegment):
    def __init__(self, id, chrom, strand, gene_start, gene_end, name='', biotype=''):
        ChromSegment.__init__(self, id, chrom, strand, gene_start, gene_end, name, biotype)
        self.transcripts = []
        self.splices = []
        self.exons = set()

    def get_fully_spliced_transcript(self):
        fully_spliced = Transcript(self.name, self.chrom, self.strand, self.start, self.end, self)
        fully_spliced.exons = list(self.exons)
        fully_spliced.valid_exons()
        return fully_spliced

    def bed_text(self): self.get_fully_spliced_transcript().bed_text()

    def gtf_text(self):
        gtf_lines = '' #TODO: make this the initial gene gtf line if I want to match Ensembl
        for transcript in self.transcripts:
            gtf_lines += transcript.gtf_text()
        return gtf_lines


class Transcript(ChromSegment):
    def __init__(self, id, chrom, strand, trans_start, trans_end, gene, name='', biotype='', start_codon_start=None, stop_codon_start=None):
        ChromSegment.__init__(self, id, chrom, strand, trans_start, trans_end, name, biotype)
        self.gene = gene
        self.splices = []
        self.exons = []
        self.start_codon_start = start_codon_start
        self.stop_codon_start = stop_codon_start

    def add_exon(self, exon):
        self.exons.append(exon)
        if True not in [exon.equals(gene_exon) for gene_exon in self.gene.exons]:
            self.gene.exons.add(exon)

    def valid_exons(self):
        exons_by_start = sorted(self.exons, key=lambda x: x.start)
        exons_by_end = sorted(self.exons, key=lambda x: x.start)

        if exons_by_start != exons_by_end: return False
        for exon_list in [exons_by_start, exons_by_end]:
            for i, exon in enumerate(exon_list):
                if i != 0 and exon.start <= exon_list[i-1].end: return False
                if i != len(exons_by_start) - 1 and exon.end >= exon_list[i+1].start: return False

        self.exons = exons_by_start
        return True

    def get_seq(self): return ''.join([str(exon) for exon in self.exons])

    def get_seq(self, start, end): return ''.join([exon.get_seq(start, end) for exon in self.exons])

    def translate(self):
        if not self.start_codon_start or not self.stop_codon_start: return self.three_frame_translate()
        seq = self.get_seq(self.start_codon_start, self.stop_codon_start)
        seq_to_translate = seq[:len(seq) - len(seq) % 3]
        return translate(seq_to_translate) if seq_to_translate and len(seq_to_translate) >= 3 else None

    def three_frame_translate(self, protein_id=''):
        frames = 3
        seq = self.get_seq()
        translations = [None for i in frames]
        for i in range(frames):
            stop_translation = len(seq) - ((len(seq) - i) % 3)
            translation = translate(seq[i:stop_translation])
            if translation:
                translations[i] = AminoAcidSequence(protein_id if protein_id else self.name + '_' + str(i), translation, self)
        return translations

    def bed_text(self):
        ChromSegment.bed_text(self, 0, len(self.exons), [len(exon) for exon in self.exons], [exon.start for exon in self.exons])

    def gtf_text(self, source='ensembl'):
        gtf_lines = ''
        ATTRIB_FIELDS = ['gene_id', 'transcript_id', 'exon_number', 'gene_name', 'gene_source', 'gene_biotype', 'transcript_name', 'transcript_source', 'transcript_biotype', 'exon_id', ]
        for i, exon in enumerate(self.exons):
            attributes = [self.gene.id, self.id, str(i), self.gene.name, source, self.gene.biotype,  self.name, source, self.biotype, exon.id]
            attributes = '; '.join([ATTRIB_FIELDS[j] + ' "' + attributes[j] + '"' for j in len(ATTRIB_FIELDS)])
            gtf_lines += exon.gtf_text() + '\t' + attributes + '\n'
        return gtf_lines


class DiExon(Transcript):
    def __init__(self, id, chrom, strand, exon1, exon2, name='', biotype=''):
        first_exon = exon1 if exon1.is_before(exon2) else exon2
        second_exon = exon2 if exon2.is_after(exon1) else exon1
        Transcript.__init__(self, id, chrom, strand, first_exon.start, second_exon.end, name, biotype) #TODO: consider the case of a chimeric junction
        self.splice = SpliceJunction(id, None, first_exon, second_exon, name, biotype)
        self.exons = [first_exon, second_exon]
        self.valid_exons()
        self.splices = [self.splice]

    def get_tryptic_splice_peptides(self):
        tryp_splice_peps = []
        translations = self.three_frame_translate()
        for i, translation in enumerate(translations):
            first_exon_stop = translation.find('*', len(self.exons[0]) / 3)
            if first_exon_stop >= 0 or not translation: continue
            second_exon_stop = translation.find('*', len(self.exons[0]) / 3)
            tstart = 0
            tend = second_exon_stop if second_exon_stop >= 0 else len(translation)
            translation.trim(tstart, tend)

            tryptic_peptides = translation.tryptic_digestion()
            for j, tryp_pep in tryptic_peptides:
                crosses_splice = tryp_pep.trans_start < self.splice.start and tryp_pep.trans_start > self.splice.end
                if crosses_splice:
                    tryp_splice_peps.append(tryp_pep)
        return tryp_splice_peps


class Exon(ChromSegment):
    def __init__(self, id, chrom, strand, exon_start, exon_end, transcript, name='', biotype=''):
        ChromSegment.__init__(self, id, chrom, strand, exon_start, exon_end, name, biotype)
        self.transcript = transcript
        self.biotype = '' #TODO: store this information from the GTF

    def get_seq(self): return str(self)

    def get_seq(self, start, end):
        if start <= self.start:
            if end >= self.end: return str(self)
            elif end >= self.start: return str(self)[:end - self.end + 1]
        elif start <= self.end:
            if end >= self.end: return str(self)[start - self.start:]
            elif end >= self.start: return str(self)[start - self.start: end - self.end + 1]
        return ''

    def get_reverse_complement(self): return reverse_complement(self.get_seq())

    def get_reverse_complement(self, start, end): return reverse_complement(self.get_seq(start, end))

    def bed_text(self): ChromSegment.bed_text(self, 0, 1, len(self), 0)

    def gtf_text(self): return '\t'.join([str(self.chrom.name)], self.biotype, 'exon', self.start, self.end, '.', self.strand, '.') #partial text, leaving off attributes


class SpliceJunction(ChromSegment):
    def __init__(self, id, chrom, strand, intron_start, intron_end, name='', biotype='', prev_exon=None, post_exon=None):
        ChromSegment.__init__(self, id, chrom, strand, intron_start, intron_end, name, biotype)
        self.prev_exon = prev_exon
        self.post_exon = post_exon

    def bed_text(self): return DiExon(self.name, self.prev_exon, self.post_exon).bed_text() if self.prev_exon and self.post_exon else ''

    def gtf_text(self): return DiExon(self.name, self.prev_exon, self.post_exon).gtf_text() if self.prev_exon and self.post_exon else ''

