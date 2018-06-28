mkdir -p db
cd db

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
md5sum taxdump.tar.gz.md5
tar -xf taxdump.tar.gz
rm taxdump.tar.gz

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz.md5
md5sum nucl_gb.accession2taxid.gz.md5
gunzip nucl_gb.accession2taxid.gz

# download everything we want from refseq -- maybe all of it!?
# pull >header lines with accessions out into a simple list
# so that we can use it to filter the acc2tax list
awk 'FNR==NR{a[$1]="y"; next}{if(a[$2]=="y"){print $0}}' refseq_accessions.txt nucl_gb.accession2taxid > nucl_gb.accession2taxid.filtered &
