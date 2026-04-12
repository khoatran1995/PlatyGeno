import os
import gzip
import requests
from tqdm import tqdm

def download_sra_subset(accession, read_limit=20000, output_dir="data"):
    """
    Downloads a small subset of reads from ENA (European Nucleotide Archive).
    This avoids downloading multi-gigabyte SRA files.
    """
    # Use ENA API to get the correct download URL
    api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=fastq_ftp"
    
    try:
        api_res = requests.get(api_url, timeout=20)
        api_res.raise_for_status()
        
        # Parse API response (TSV format)
        lines = api_res.text.strip().split("\n")
        if len(lines) < 2:
            print(f"Error: ENA API returned no data for {accession}")
            return
            
        # Get the fastq_ftp field (index 1)
        data_line = lines[1].split("\t")
        if len(data_line) < 2:
             print(f"Error: ENA API response format unknown for {accession}")
             return
             
        ftp_links = data_line[1].split(";")
        # Take the first file (often _1.fastq.gz for paired-end)
        ftp_path = ftp_links[0]
        
        # Translate FTP to HTTP for better 'requests' compatibility
        url = f"http://{ftp_path}"
        
    except Exception as e:
        print(f"Error fetching metadata from ENA: {e}")
        # Fallback to hardcoded pattern if API fails
        prefix = accession[:6]
        subdir = accession[-3:]
        url = f"http://ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}/{subdir}/{accession}/{accession}_1.fastq.gz"
    
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, f"{accession}_subset.fastq")
    
    print(f"Requesting stream from ENA: {url}")
    print(f"Target: {read_limit} reads...")

    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        # Decompress the stream on the fly
        # We read chunks and use a Gzip decompressor
        import zlib
        decompressor = zlib.decompressobj(32 + zlib.MAX_WBITS) # handles gzip headers
        
        reads_count = 0
        buffer = ""
        
        with open(out_file, "w") as f:
            for chunk in response.iter_content(chunk_size=1024 * 1024): # 1MB chunks
                if not chunk:
                    break
                
                try:
                    data = decompressor.decompress(chunk).decode("utf-8")
                    buffer += data
                except Exception:
                    # If we hit a decompression error, it might be the end of the stream
                    break
                    
                lines = buffer.split("\n")
                # Leave the last (potentially incomplete) line in the buffer
                buffer = lines.pop()
                
                for i in range(0, len(lines), 4):
                    if i + 3 < len(lines):
                        # Write the 4-line FASTQ record
                        f.write("\n".join(lines[i:i+4]) + "\n")
                        reads_count += 1
                        
                        if reads_count >= read_limit:
                            break
                
                if reads_count >= read_limit:
                    break
                    
        print(f"Downloaded {reads_count} reads to {out_file}")
        
    except Exception as e:
        print(f"Failed to download subset: {e}")
        print("Tip: Ensure you have a stable internet connection or try a different accession.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Download a small subset of an exotic metagenome.")
    parser.add_argument("accession", type=str, help="SRR Accession (e.g., SRR23196177)")
    parser.add_argument("--limit", type=int, default=20000, help="Number of reads to download")
    parser.add_argument("--out", type=str, default="data", help="Output directory")
    
    args = parser.parse_args()
    download_sra_subset(args.accession, args.limit, args.out)
