import os
import argparse
from ftplib import FTP

def ftp_download_sra(accession, folder="data"):
    server = "ftp.sra.ebi.ac.uk"
    
    # ENA FTP directory structure: /vol1/fastq/SRR546/029/SRR5462529/
    prefix = accession[:6]
    suffix = f"0{accession[-2:]}" if len(accession) > 9 else f"00{accession[-1:]}"
    # Note: Modern ENA structure can vary, we'll try to find it
    path = f"/vol1/fastq/{prefix}/{suffix}/{accession}/"
    
    os.makedirs(folder, exist_ok=True)
    
    print(f"Connecting to ENA FTP: {server}")
    try:
        ftp = FTP(server)
        ftp.login() # Anonymous login
        
        print(f"Navigating to: {path}")
        ftp.cwd(path)
        
        # List files to be sure
        remote_files = ftp.nlst()
        print(f"Found files: {remote_files}")
        
        for filename in remote_files:
            if not filename.endswith(".fastq.gz"):
                continue
                
            output_path = os.path.join(folder, filename)
            print(f"Downloading {filename}...")
            
            with open(output_path, 'wb') as f:
                ftp.retrbinary(f"RETR {filename}", f.write)
            
            print(f"Download Complete: {filename}")
            
        ftp.quit()
        
    except Exception as e:
        print(f"FTP Error: {e}")
        print("Tip: If directory not found, the accession might be in a different volume. Check ENA portal.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download SRA data via FTP.")
    parser.add_argument("accession", type=str, help="SRR Accession")
    parser.add_argument("--out", type=str, default="data", help="Output directory")
    
    args = parser.parse_args()
    ftp_download_sra(args.accession, args.out)
