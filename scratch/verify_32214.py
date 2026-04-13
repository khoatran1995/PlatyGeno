from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time

seq_60 = "TATCTTTTGGCGCATATTCTCTGCCCCGTTACAATCTGCGTTCGCTACCAATTCGCACGA"
seq_108 = "TATCTTTTGGCGCATATTCTCTGCCCCGTTACAATCTGCGTTCGCTACCAATTCGCACGACGAACACACGTACAAACCACGATGCTTTCGGTTCGATTTCGTCGTATCAC"

def check_blast(seq, label):
    print(f"Checking {label} ({len(seq)} bp)...")
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", seq, format_type="XML")
        blast_record = NCBIXML.read(result_handle)
        if len(blast_record.alignments) > 0:
            top_hit = blast_record.alignments[0]
            hsp = top_hit.hsps[0]
            identity = (hsp.identities / hsp.align_length) * 100
            print(f"  Hit: {top_hit.title[:80]}")
            print(f"  Identity: {identity:.2f}%")
            print(f"  E-value: {hsp.expect}")
        else:
            print("  No hits found in 'nt' database.")
    except Exception as e:
        print(f"  Error: {e}")

check_blast(seq_60, "Precision Snippet")
time.sleep(5)
check_blast(seq_108, "Consensus Assembly")
