# module name: datamanagers
# main program: samplespecificdbgenerator

__author__ = "anthony"
__date__ = "$Oct 27, 2015 2:54:44 PM$"

import genemodel
import variantcalls
from lxml import etree as et

#Expect to make a new database if no path is given. Otherwise, manage the manimpulation of a current one.
class ProteinXmlManager:
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
            headerInfo = [x for x in [accession, name, chromosome, geneId, transcriptId, addedInfo, score] if
                          x != None]  # remove None elements
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
class FastaManager:
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