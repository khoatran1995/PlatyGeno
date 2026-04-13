import urllib.request
import os

def download_full_sra(accession, folder="data/full_benchmarking"):
    # ENA FTP URLs for _1 and _2
    urls = [
        f"https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/077/{accession}/{accession}_1.fastq.gz",
        f"https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/077/{accession}/{accession}_2.fastq.gz"
    ]
    
    os.makedirs(folder, exist_ok=True)
    
    for url in urls:
        filename = os.path.basename(url)
        output_path = os.path.join(folder, filename)
        
        print(f"Downloading FULL SRA: {url}")
        try:
            # Note: reporthook only shows progress for the current file
            urllib.request.urlretrieve(url, output_path, reporthook=progress_report)
            print(f"\nDownload Complete: {filename}")
        except Exception as e:
            print(f"\nError during download {filename}: {e}")

def progress_report(block_num, block_size, total_size):
    downloaded = block_num * block_size
    progress = (downloaded / total_size) * 100 if total_size > 0 else 0
    if block_num % 1000 == 0:
        # Using a more robust printing for Windows
        print(f"   Downloaded: {downloaded / 1024 / 1024:.1f} MB ({progress:.1f}%)", end='\r')

if __name__ == "__main__":
    download_full_sra("SRR23196177")
