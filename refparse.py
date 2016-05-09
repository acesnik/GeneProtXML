# module name: datamanagers
# main program: samplespecificdbgenerator

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:54:44 PM$"

import re
import genemodel
import numpy as np
from lxml import etree as et


class FileManager:
    def __init__(self, file=None):
        pass

#Expect to make a new database if no path is given. Otherwise, manage the manimpulation of a current one.
class ProteinXmlManager(FileManager):
    def __init__(self, xml=None):
        HTML_NS = "http://uniprot.org/uniprot"
        XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
        NAMESPACE_MAP = {None: HTML_NS, "xsi": XSI_NS}
        self.UP = '{' + HTML_NS + '}'

        self.db, self.root = None, None
        if not xml:
            self.root = et.Element(self.UP + 'uniprot', nsmap=NAMESPACE_MAP)
            self.db = et.ElementTree(self.root)
        elif isinstance(xml, et._ElementTree): self.db = xml
        elif isinstance(xml, et._Element): 
            print 'Error: use the ElementTree in ProteinXmlManager, not an Element'
            exit(2)
        else: self.db = et.parse(xml)
        self.root = self.db.getroot()

    def lxml_entry(self, object):
        if isinstance(object, genemodel.SequenceVariant):
            pass
        if isinstance(object, genemodel.AminoAcidSequence):
            entry = self.blank_xml_entry()
            return entry
        if isinstance(object, genemodel.Indel):
            pass
        if isinstance(object, genemodel.Modification):
            pass

    def condense_xml_entry(self, entry):
        for element in entry:
            if element.tag not in [self.UP + 'protein', self.UP + 'accession', self.UP + 'name', self.UP + 'gene',
                                   self.UP + 'organism', self.UP + 'proteinExistence', self.UP + 'depth',
                                   self.UP + 'sequence', self.UP + 'feature', self.UP + 'dbReference']:
                entry.remove(element)
            elif element.get('type') != 'Ensembl' and element.tag == UP + 'dbReference':
                entry.remove(element)
            elif element.tag == self.UP + 'organism':
                for field in element:
                    if field.tag != self.UP + 'name':
                        element.remove(field)
            elif element.tag == self.UP + 'protein':
                for name in element:
                    if name.tag != self.UP + 'recommendedName':
                        element.remove(name)

    def condense_xml_entries(self, root):
        for entry in root:
            self.condense_xml_entry(entry)

    def blank_xml_entry(self):
        entry = et.Element(self.UP + 'entry', dataset="Ensembl")
        accession = et.SubElement(entry, self.UP + 'accession')
        name = et.SubElement(entry, self.UP + 'name')
        fullName = et.SubElement(et.SubElement(et.SubElement(entry, self.UP + 'protein'), self.UP + 'recommendedName'), self.UP + 'fullName')
        gene = et.SubElement(entry, self.UP + 'gene')
        geneName1 = et.SubElement(gene, self.UP + 'name', type="coords")
        organism = et.SubElement(entry, self.UP + 'organism')
        organismName1 = et.SubElement(organism, self.UP + 'name', type="scientific")
        organismName2 = et.SubElement(organism, self.UP + 'name', type="common")
        organismName1.text = "Homo sapiens"
        organismName2.text = "Human"
        proteinExist = et.SubElement(entry, self.UP + 'proteinExistence', type="evidence at transcript level")
        sequence = et.SubElement(entry, self.UP + 'sequence', version="1", fragment="single")
        return entry

    def set_sequence(self, entry, sequence):
        entry.find(self.UP + 'sequence').text = sequence.replace('\n', '').replace('\r', '')

    def get_all_sequence_elements(self):
        return self.root.findall('.//'+self.UP+'sequence')

    def get_all_sequences(self):
        return [e.text for e in self.get_all_sequence_elements()]

    # Returns information in an XML entry in fasta duplet: (header, sequence)
    def xml_to_fasta(self, entry):
        header = ">"
        if entry.tag == self.UP + 'copyright': return None
        if entry.get('dataset') == 'Ensembl':
            accession, name, geneInfo, chromosome, geneId, transcriptId, addedInfo = None, None, None, None, None, None, None
            accession = entry.find(self.UP + 'accession').text
            name = entry.find(self.UP + 'name').text
            geneInfo = entry.find(self.UP + 'gene')
            if geneInfo != None:
                geneInfo.getiterator(self.UP + 'name')
                for item in geneInfo:
                    if item.get('type') == 'coords': chromosome = item.text
                    if item.get('type') == 'primary':
                        geneId = 'gene:' + item.text
                        transcriptId = 'transcript:' + item.find(self.UP + 'transcript').text
            addedInfo = entry.find(self.UP + 'protein').find(self.UP + 'recommendedName').find(self.UP + 'fullName').text
            score = entry.find(self.UP + 'proteinExistence').find(self.UP + 'depth')
            score = 'depth:' + score.get('reads') if score != None else ''
            headerInfo = [x for x in [accession, name, chromosome, geneId, transcriptId, addedInfo, score] if x != None]  # remove None elements
            header += ' '.join(headerInfo)
        else:
            if entry.get('dataset') == 'Swiss-Prot': database = 'sp'
            else: database = 'tr'
            accession = entry.find(self.UP + 'accession').text
            name = entry.find(self.UP + 'name').text
            accession = '|'.join([database, accession, name])
            organism = entry.find(self.UP + 'organism')
            if organism != None: organism = 'OS=' + organism.find(self.UP + 'name').text
            geneName = entry.find(self.UP + 'gene')
            if geneName != None: geneName = 'GN=' + geneName.find(self.UP + 'name').text
            headerInfo = [x for x in [accession, organism, geneName] if x != None]
            header += ' '.join(headerInfo)
        return [header.replace('\n', '').replace('\r', ''), entry.find(self.UP + 'sequence').text.replace('\n', '').replace('\r', '')]
    
    def write_to_fasta(self, fasta_out):
        fasta_out.write('\n'.join(['\n'.join(self.xml_to_fasta(entry)) for entry in self.root]))
            

#Currently have no need to create fasta objects. Only takes in a fasta and manages it.
#Make flexible, for both genome and protein at least
class FastaManager(FileManager):
    def __init(self, fasta=None):
        self.fasta = fasta

    # Get header and seq from protein fasta using transcript or protein ensembl accession
    def read_fasta_to_xml(self, refFasta, db=None):
        xml_manager = ProteinXmlManager(db)
        seq = ""
        line = refFasta.readline().strip()
        while line != "":
            if line.startswith(">"):
                line = line.split()
                acc, seqtype, chromosome = line[0][1:], line[1], line[2]
                geneId, transcriptId, addedInfo = line[3].split(':')[1], line[4].split(':')[1], ' '.join(line[5:])
                line = refFasta.readline().strip()
                while not line.startswith(">"):
                    if line == "": break
                    seq += line
                    line = refFasta.readline().strip()
                if seq.find('*') < 0: xml_manager.enter_seqvar(root, acc, seqtype, chromosome, addedInfo, '', '',geneId, transcriptId, seq) #TODO: implement this
                seq = ""

    # Reads the headers and sequences of a fasta file into RAM
    def read_fasta(self, fasta):
        header_seq_pair = ([], [])  # headers, #sequences
        line = fasta.readline()
        while line.startswith('#'):
            line = fasta.readline()
        sequence = ""
        while line != "":
            if line.startswith(">"):
                header_seq_pair[0].append(line)
                line = fasta.readline()
                while not line.startswith(">"):
                    if line == "": break
                    sequence += line
                    line = fasta.readline()
                    header_seq_pair[1].append(sequence.replace('\n', '').replace('\r', ''))
                sequence = ""
        fasta.close()
        return header_seq_pair


#Takes in a GTF file and places it into the genemodel objects
class GtfManager(FileManager):
    def __init__(self, gtf=None):
        self.gtf = gtf


#Takes in variant call information and places it into the genemodel objects
class VcfManager(FileManager):
    def __init__(self, vcf_path, chromosomes):
        self.vcf_path = vcf_path
        self.chromosomes = chromosomes
        self.header = []
        self.field_header = []
        self.individuals = []
        self.haplotype_stack = [] #Used to handle haplotype grouping when parsing variant call lines

    def read_vcf(self):
        for line in open(self.vcf_path, 'r'):
            if line.startswith('##'):
                self.header.append(line)
            elif line.startswith('#CHROM'):
                self.field_header = line[1:].split('\t')
                individual_names = self.field_header[10:] if len(self.field_header) >= 10 else ['individual1']
                self.individuals = [genemodel.Individual(name, self.chromosomes) for name in individual_names]
                self.haplotype_stack = [[] for name in individual_names]
            else:
                self.parse_vcf_line(line)

    # This method considers the phasing information recorded in the GT field for each sample to initialize haplotypes
    # with the start and stop loci indicated by the variant calls.
    # Example: chromosome 1, position 1 through 8 have the genotypes [0/1, 1|1, 0|1, 0/1, 0|1, 0/1, 0|1, 0|1]
    # Note: 0/1 or 1/0 are heterozygous alleles unphased with the previous variant call
    # Note: 0|1 or 1|0 are heterozygous alleles phased with the previous variant call
    # Note: 0|0 or 1|1 are homozygous alleles and are always phased with the above
    # Note: Alternate alleles are indicated by a genotype > 0. If there is a list of alternate alleles,
    # they are specified by 1, 2, and so on.
    # Note: Haploid chromosomes (MT, Y, and X in males) should be declared unphased, and I plan to ignore "heterozygous"
    # haploid variant calls.
    def parse_variant(self, line):
        fields = line.split('\t')
        (chrom, position, id, reference, alternates, qual, filter, info, format) = [fields[i] if i < len(fields) else None for i in range(9)]
        individual_info_list = fields[10:] if len(fields) >= 10 else None
        alternates = alternates.split(',') #There can be multiple alternates
        qual = float(qual)
        depth = -1
        allele_frequency = -1
        genotype = None
        if info:
            for info_item in info.split(';'):
                if info_item.find('=') < 0: return
                (key, val) = info_item.split('=', 1)
                if key.startswith('AF'): allele_frequency = float(val)
                if key.startswith('DP'): depth = int(val)

        seqvars = [genemodel.SequenceVariant(chrom, position, id, reference, alternate, qual, allele_frequency, depth) for alternate in alternates]

        # No individual-specific information. Just append the seqvar(s) to a generic individual
        if not format: self.individuals[0].sequence_variants += seqvars

        # Individual-specific information available
        if format:
            for i, individual_info in enumerate(individual_info_list):
                self.individuals[i].sequence_variants += seqvars
                for indiv_info_item in individual_info.split(':'):
                    for f, format_field in enumerate(format.split(':')):
                        if format_field.startswith('DP'): depth = int(indiv_info_item[f]) #overwrite info standard field, since this one includes filtering criteria
                        if format_field.startswith('GT'): genotype = indiv_info_item[f]
                genotype = list(genotype)
                ploid1_allele_ref, phase, ploid2_allele_ref = int(genotype[0]), genotype[1], int(genotype[2]) # 0 refers to reference. >0 refers to one of the alternates
                ploid1_allele = seqvars[int(ploid1_allele) - 1] if ploid1_allele else None
                ploid2_allele = seqvars[int(ploid2_allele) - 1] if ploid2_allele else None
                is_phased = phase == '|'

                current_haplotype = self.individuals[i].local_haplotypes[-1]
                if not is_phased:
                    self.close_haplotype(current_haplotype, i)
                    self.individuals[i].local_haplotypes.append(genemodel.LocalHaplotype(chrom, position))
                self.haplotype_stack[i].append(ploid1_allele, ploid2_allele)
                current_haplotype.update_end(max(position, ploid1_allele.end if ploid1_allele else 0, ploid2_allele.end if ploid2_allele else 0))

    def close_haplotype(self, haplotype, individual_index):
        for ploid1_seqvar, ploid2_seqvar in self.haplotype_stack[individual_index]:
            haplotype.add(ploid1_seqvar, ploid2_seqvar)


#Takes in a bed file, such as a Tophat splice junction bed file
class BedManager(FileManager):
    def __init__(self, bed=None):
        self.bed = bed


#Takes in annotations from slncky, such as the filtered_info file
class SlnckyManager(FileManager):
    def __init__(self):
        pass


#Takes in a STAR SJ.out.tab file specifying splice junctions
class StarManager(FileManager):
    def __init__(self):
        pass


class ProBamManager(FileManager):
    def __init__(self):
        pass


def add_unified(root, newId, rootIndex, acc, seqtype, chromosome, biotypes, geneId, transcriptId, seq):
    if newId:
        entry = et.SubElement(root, UP+'entry', dataset="Ensembl")
        accession = et.SubElement(entry, UP+'accession')
        accession.text = str(newId)
    else: entry = root[rootIndex]

    accession = et.Element(UP+'accession')
    accession.text = acc
    a = entry.findall(UP+'accession')
    if len(a) > 0: a[-1].addnext(et.Element(UP+'accession'))
    else: entry.append(accession)

    name = et.Element(UP+'name')
    name.text = seqtype
    n = entry.findall(UP+'name')
    if len(n) > 0: n[-1].addnext(et.Element(UP+'name'))
    else: entry.append(name)

    if newId:
        fullName = et.SubElement(et.SubElement(et.SubElement(entry, UP+'protein'), UP+'recommendedName'), UP+'fullName')
        fullName.text = 'referenceProtein' + str(newId)
        gene = et.SubElement(entry, UP+'gene')
    else:
        gene = et.Element(UP+'gene')
        entry.findall(UP+'gene')[-1].addnext(gene)
    geneName1 = et.SubElement(gene, UP+'name', type="coords")
    geneName1.text = chromosome
    geneName2 = et.SubElement(gene, UP+'name', type="primary")
    geneName2.text = geneId
    gene_biotype = et.SubElement(gene, UP+'gene_biotype')
    gene_biotype.text = biotypes.split(' ')[0]
    transcript = et.SubElement(gene, UP+'transcript', type="primary")
    transcript.text = transcriptId
    transcript_biotype = et.SubElement(gene, UP+'transcript_biotype')
    transcript_biotype.text = biotypes.split(' ')[1]

    if newId:
        organism = et.SubElement(entry, UP+'organism')
        organismName1 = et.SubElement(organism, UP+'name', type="scientific")
        organismName2 = et.SubElement(organism, UP+'name', type="common")
        organismName1.text = "Homo sapiens"
        organismName2.text = "Human"
        proteinExist = et.SubElement(entry, UP+'proteinExistence', type="evidence at transcript level")
        sequence = et.SubElement(entry, UP+'sequence', version="1", fragment="single")
        sequence.text = seq

#TODO: oh that's right. There are a bunch of duplicate sequences in the Ensembl database
def ensembl_entry(uniqSeqs, root, acc, seqtype, chromosome, biotypes, geneId, transcriptId, seq):
    found = False
    for i, s in enumerate(uniqSeqs):
        if seq == s:
            found = True
            add_unified(root, '', i, acc, seqtype, chromosome, biotypes, geneId, transcriptId, seq)
            break
    if not found:
        uniqSeqs.append(seq)
        add_unified(root, str(len(uniqSeqs)), -1, acc, seqtype, chromosome, biotypes, geneId, transcriptId, seq)