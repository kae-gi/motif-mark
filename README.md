# motif-mark

### *A program that visualizes to-scale motifs on sequences.* 

Given a FASTA file and motif file (file that contains the desired motifs), this program will illustrate the gene and its introns, exons, and location of motifs. Inputs from FASTA must in form of `<intron><EXON><intron>` sequences. Output is of type .png. 
It can handle:
- Input FASTA file (seqs ≤1000 bases) 
- Input motifs file (≤10 bases each, one motif per line in a text file)
- Motifs with ambiguous nucleotides
- Multiple sequences
- Multiple motifs (currently 7 max color wise, but technically can handle more with a larger color palette) 

### Example:
Below is an example figure that utilizes the provided example FASTA and motif file. The key shows the motif color code. The thicker black box is the exon, and the thinner black line is the intron. 

![Figure_1](https://user-images.githubusercontent.com/81830809/223919799-060ed3e4-c599-43c1-bc99-45d3a6bcdcc0.png)

### Usage and Requirements
- **Usage:** `./motif-mark-oop.py -f <fasta file> -m <motifs file>`
- [Pycairo package](https://pycairo.readthedocs.io/en/latest/) must be installed to run ("It depends on cairo >= 1.15.10 and works with Python 3.7+...")
- This program currently runs on Python v3.10.8

