from time import sleep
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from bioservices import UniProt
import pandas as pd
import io
from from_SGD import add_sgd_data


class GetInfo:
    '''
    From PAS_... (Pichia pastoris gene symbol) get:
        - gene Entrez ID
        - gene name
        - protein accession number (XP...)
        - protein function
    '''

    EMAIL = "t.ianshina99@gmail.com"

    def __init__(self, pichia_genes: list):
        self.pichia_genes = pichia_genes
        self.result = {}

    def get_info_from_pas(self):
        if len(self.pichia_genes) == 0:
            raise ValueError("Gene list is empty")

        Entrez.email = self.EMAIL

        for i, gene in enumerate(self.pichia_genes):
            print(
                f"{i + 1}/{len(self.pichia_genes)} - Process xp_from_pas for gene - '{gene}'"
            )

            if self.result.get(gene) is not None:
                raise KeyError(f"Такой ген ({gene}) уже есть")

            self.result[gene] = []

            with Entrez.esearch(db="gene", term=gene) as handle:
                search_record = Entrez.read(handle)
            sleep(0.2)
            gene_id = search_record["IdList"][0]

            with Entrez.efetch(db="gene", id=gene_id, retmode="xml") as handle:
                gene_record = Entrez.read(handle)

            protein_func = gene_record[0]["Entrezgene_prot"]["Prot-ref"]['Prot-ref_desc']

            protein_info = gene_record[0]["Entrezgene_locus"][0][
                "Gene-commentary_products"
            ][0]["Gene-commentary_products"][0]
            assert protein_info["Gene-commentary_type"].attributes["value"] == "peptide"
            prot_accession = protein_info["Gene-commentary_accession"]

            # uniprot
            uni = UniProt()
            query = gene

            cols = "gene_names, protein_name, go_p, go_f, go_c"
            res_uniprot = uni.search(query, frmt="tsv", columns=cols)
            sleep(0.2)
            df = pd.read_table(io.StringIO(res_uniprot), header=0)

            df.to_csv('file.csv')
            if res_uniprot is None:
                gene_name_short = None
                gene_name_full = None
                go_p = None
                go_f = None
                go_c = None
            pd.options.display.max_colwidth = 1000

            gene_name_short = df['Gene Names'].to_string(index=False).split()[0]
            if gene_name_short == gene:
                gene_name_short = None
            gene_name_full = df['Protein names'].to_string(index=False).split('\n')[0].strip()

            for col in (gene_id, gene_name_short, gene_name_full,
                        prot_accession, protein_func):
                self.result[gene].append(col)

            for g in ('biological process', 'molecular function', 'cellular component'):
                go = df[f'Gene Ontology ({g})'].to_string(index=False).replace('NaN', '').replace(';', '\n').strip()
                strings = pd.Series((i.strip() for i in go.split('\n')))
                go = strings.drop_duplicates().to_string(index=False).strip()
                print(go)
                self.result[gene].append(go)
                sleep(1)

            yield gene, prot_accession


class Blast(GetInfo):
    def protein_blast(self):
        '''
        Installed BLAST+ is necessary. The path to blastp.exe file must be written in 'blastp' variable.
        Database files must be located in the current directory. Create database files with command line:
            makeblastdb -in {protein_sequence.fasta} -parse_seqids -dbtype prot -out {database_name}
        '''

        blastp = r"C:\Program Files\NCBI\blast-2.14.0+\bin\blastp.exe"

        for gene, prot in self.get_info_from_pas():
            print(f"\t     -     Blast for gene     - '{gene}'")

            # find fasta sequence
            handle = Entrez.efetch(db='protein', id=prot, rettype='fasta', retmode='text')
            sleep(0.2)
            record = SeqIO.read(handle, 'fasta')
            seq = str(record.seq)

            # create fasta file
            my_seq = SeqRecord(Seq(seq), id=prot)
            with open("seq.fasta", "w") as output_handle:
                SeqIO.write(my_seq, output_handle, "fasta")

            # BLAST
            cline = NcbiblastpCommandline(blastp, query="seq.fasta", db="S288Cdb", outfmt=5,
                                          evalue=0.001, out="results.xml")
            stdout, stderr = cline()
            sleep(1)
            self.parse_blast_xml(gene)
        return self.result

    def parse_blast_xml(self, gene):
        """
        Takes selected data from blast.xml file:
                SC_Accession,
                SC_description,
                Query_Cover,
                E_value,
                Per_Ident,
                Acc_Len
            )
        """

        for record in NCBIXML.parse(open("results.xml")):
            if len(record.alignments) == 0:
                sc_accession = None
                sc_description = None
                accession_len = None
                query_cover = None
                e_value = None
                percent_identity = None

            else:
                filtered_alignments = list(
                    filter(lambda x: "NP" in x.accession, (i for i in record.alignments))
                )
                if len(filtered_alignments) == 0:
                    align = record.alignments[0]
                else:
                    align = filtered_alignments[0]

                sc_accession = align.accession
                sc_description = align.hit_def
                accession_len = align.length

                for hsp in align.hsps:
                    percent_identity = hsp.identities / hsp.align_length * 100
                    query_cover = (hsp.query_end - hsp.query_start + 1) / record.query_length * 100
                    e_value = hsp.expect

                sgd = add_sgd_data(sc_description)

        for sc_col in (sc_accession, sc_description, sgd, query_cover,
                       e_value, percent_identity, accession_len):
            self.result[gene].append(sc_col)
