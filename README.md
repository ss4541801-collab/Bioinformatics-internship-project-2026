# 🧬 BRCA1 Gene Analysis Pipeline
### DNA Sequence Analysis & Open Reading Frame Detection using Python & Biopython

> **CBIRT Summer Internship Project — May to June 2026**  
> Parul Institute of Technology, Parul University, Vadodara

---

## 📌 What is this project?

This is a Python-based bioinformatics pipeline that automatically fetches the **BRCA1 gene** (Breast Cancer Gene 1) from the NCBI database and performs a full sequence analysis — from raw nucleotides all the way to a translated protein sequence.

BRCA1 is one of the most studied human genes. Mutations in it can increase the risk of breast cancer up to **72%** and ovarian cancer up to **44%**. This project explores the gene computationally using free, open-source tools.

---

## ⚙️ What the pipeline does

| Step | Task | Output |
|------|------|--------|
| 1 | Fetch BRCA1 mRNA sequence from NCBI | `DNA.fasta` |
| 2 | Calculate GC content | Printed to console |
| 3 | Count and visualise nucleotides (A, T, G, C) | `nucleotide_composition.png` |
| 4 | Find all ORFs across 6 reading frames | `Forward_ORF.txt`, `Backward_ORF.txt` |
| 5 | Translate the longest ORF into protein | `Translation_seq.txt` |
| 6 | Generate reverse complement sequence | `Reverse_compl_seq.txt` |

---

## 📊 Results (NCBI Accession: NM_007294)

| Analysis | Result |
|----------|--------|
| Sequence Length | 7,088 nucleotides |
| GC Content | 41.77% ✅ matches human genome average |
| Nucleotide Counts | A: 2116 &nbsp; T: 2074 &nbsp; G: 1474 &nbsp; C: 1424 |
| Total ORFs Found | 98 (across all 6 reading frames) |
| Longest ORF | 5,592 nucleotides (forward strand, frame +2) |
| Protein Length | **1,863 amino acids** — exact match to known BRCA1 protein ✅ |

---

## 🛠️ Tools & Libraries Used

- **Python 3**
- **Biopython** — `Entrez` for NCBI fetch, `SeqIO` for FASTA, `Seq` for analysis
- **Matplotlib** — nucleotide composition bar chart
- **NCBI RefSeq** — data source (NM_007294)

---

## 🚀 How to run this project

**1. Clone the repository**
```bash
git clone https://github.com/ss4541801-collab/Bioinformatics-internship-project-2026.git
cd Bioinformatics-internship-project-2026
```

**2. Install dependencies**
```bash
pip install biopython matplotlib
```

**3. Run the script**
```bash
python main.py
```

**4. Enter the accession number when prompted**
```
Enter the accession number: NM_007294
```

> 💡 You can enter **any NCBI accession number** — the pipeline works for any gene, not just BRCA1!

---

## 📁 Output Files

After running, these files will appear in your folder:

| File | Description |
|------|-------------|
| `DNA.fasta` | Raw BRCA1 mRNA sequence in FASTA format |
| `nucleotide_composition.png` | Bar chart of A, T, G, C counts |
| `Forward_ORF.txt` | All ORFs found on the forward strand |
| `Backward_ORF.txt` | All ORFs found on the reverse complement strand |
| `Translation_seq.txt` | Protein sequence of the longest ORF |
| `Reverse_compl_seq.txt` | Full reverse complement sequence |

---

## 🧠 Key Concepts Used

- **GC Content** — percentage of G and C bases; indicates sequence stability
- **Open Reading Frame (ORF)** — region starting with ATG and ending at a stop codon (TAA, TAG, TGA)
- **6 Reading Frames** — 3 forward (+1, +2, +3) + 3 reverse complement (−1, −2, −3)
- **Translation** — converting the longest ORF into an amino acid (protein) sequence
- **Reverse Complement** — the complementary strand read in reverse (as it would be in real DNA)

---

## 🏛️ Internship Details

| | |
|--|--|
| **Organisation** | Centre of Bioinformatics Research & Technology (CBIRT) |
| **Duration** | 4 May 2026 – 30 June 2026 (2 months) |
| **Mentors** | Ms. Anchal Negi & Dr. Tamanna Anwar |
| **Focus** | Python for Bioinformatics |
| **Institute** | Parul Institute of Technology, Parul University, Vadodara |

---

## 👤 Author

**Shubham** — B.Tech Biotechnology, Parul University  
📧 ss4541801@gmail.com  
🔗 [GitHub](https://github.com/ss4541801-collab)

---

*Made with 🧬 Python, Biopython, and a lot of curiosity.*
