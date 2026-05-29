# TUGAS PROYEK IF3211 KOMPUTASI DOMAIN SPESIFIK - Kelompok 16 Julius Genetics

Proyek ini merupakan implementasi algoritma komputasi bioinformatika untuk menganalisis kemiripan sekuens genom SARS-CoV-2 dan mengklasifikasikan varian virus tersebut berdasarkan hasil alignment dengan implementasi dua algoritma sequence alignment klasik, yaitu Needleman-Wunsch untuk global alignment dan Smith-Waterman untuk local alignment. Analisis dilakukan pada sekuens gen Spike yang diekstraksi dari genom SARS-CoV-2 menggunakan pendekatan coordinate slicing berdasarkan posisi gen pada genom referensi.

Hasil dari proses perhitungan skor alignment terhadap sekuens referensi digunakan oleh program untuk melakukan klasifikasi varian SARS-CoV-2 dan membangun pohon filogenetik menggunakan Biopython untuk memvisualisasikan hubungan kekerabatan antar sekuens. Hasil akhir dari program berupa suatu klasifikasi varian pada sampel uji terhadap referensi sampel varian serta visualisasi pohon filogenetik dalam format Newick, PNG, dan SVG.

# Dataset

Dataset yang digunakan berasal dari NCBI GenBank untuk varian Sars-Cov-2 dan terdiri atas sekuens referensi untuk proses klasifikasi serta sekuens sampel uji untuk proses evaluasi.

### Sekuens Referensi

| Varian | Accession |
|---------|-----------|
| Alpha | OQ898928.1 |
| Beta | OR578390.1 |
| Delta | PV848452.1 |
| Omicron | PX937219.1 |

### Sekuens Sampel Uji

| Varian Aktual | Accession |
|--------------|-----------|
| Alpha | MZ305031.1 |
| Beta | OR353131.1 |
| Delta | PV848449.1 |
| Omicron | PP846150.1 |

## Struktur Prooyek
- `main.py` menjalankan demo end-to-end.
- `src/alignments.py` berisi implementasi alignment.
- `src/classifier.py` berisi klasifikasi varian.
- `src/phylo.py` berisi pembuatan pohon.
- `data/refs.fasta` berisi sekuens referensi per varian.
- `data/sequences.fasta` berisi sekuens uji.
- `output/tree.nwk` hasil pohon (Newick).
- `output/tree.png` dan `output/tree.svg` visualisasi pohon.

## How To Run
1. Clone Repository
```bash
git clone https://github.com/julius123123/Tubes1-KDS-KLPK-12.git
```
2. Buat dan aktifkan Virtual Environment
```bash
python -m venv .venv
.\.venv\Scripts\activate
```
3. Install Dependencies
```bash
pip install -r requirements.txt
```
3. Run the Program (di direktori utama)
```bash
python main.py
```

## Anggota Kelompok

| Nama | NIM |
|------|-----|
| Muhammad Aufa Farabi | 13523023 |
| Joel Hotlan Haris Siahaan | 13523025 |
| Julius Arthur | 13523030 |
| Ferdinand Gabe Tua Sinaga | 13523051 |

