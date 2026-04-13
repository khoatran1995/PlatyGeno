import os
import urllib.request
import gzip
import shutil

def download_subset(ftp_url, output_path, num_reads=50000):
    """Downloads a subset of reads from a gzipped FASTQ via FTP stream."""
    num_lines = num_reads * 4
    print(f"Streaming {num_reads} reads from: {ftp_url}")
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    try:
        with urllib.request.urlopen(ftp_url) as response:
            with gzip.open(response, 'rt') as gfile:
                with open(output_path, 'w') as outfile:
                    for i, line in enumerate(gfile):
                        if i >= num_lines:
                            break
                        outfile.write(line)
                        if (i + 1) % 10000 == 0:
                            print(f"   [Progress] Processed {int((i+1)/4)} reads...")
                            
        print(f"Subset complete: {output_path}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    url = "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/008/SRR1040858/SRR1040858.fastq.gz"
    out = "data/SRR1040858/bonney_50k.fastq"
    download_subset(url, out)
