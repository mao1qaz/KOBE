 <p align="center">
  <img src="docs/images/KOBE_logo.png" alt="KOBE Logo" width="240">
</p>

<h1 align="center">KOBE</h1>

<p align="center">
  <em>Kit for orthology-based editing analysis – 8 modules in 1</em>
</p>

<p align="center">
  A command-line toolkit for A-to-I RNA editing identification and analysis.
</p>

---

## 1. About KOBE

KOBE is a command-line toolkit for identifying and analyzing adenosine-to-inosine (A-to-I) RNA editing. It integrates de novo RNA editing site identification, orthology-based RNA editing analysis, orthology-guided site recovery, differential RNA editing analysis, and downstream visualization and data processing. Its modular design supports flexible analyses across samples and species.

## 2. Main Features

KOBE provides eight major functionalities:

1. **De novo RNA Editing Site Identification:** Identifies candidate RNA editing sites from RNA sequencing data, with genomic DNA sequencing data used to distinguish RNA editing events from genomic variants. Identified sites can be further annotated to determine their genomic and functional characteristics.

2. **Orthology-based RNA Editing Analysis:** Analyzes RNA editing sites across species based on orthologous gene or transcript relationships, enabling the identification and comparison of conserved and species-specific editing sites.

3. **Conservation Visualization:** Generates visual summaries of conserved RNA editing sites and their distribution across species, facilitating comparisons of conserved and species-specific editing patterns.

4. **Orthology-guided Site Recovery:** Recovers potential RNA editing sites in a target species using homologous or orthologous sequence information from other species, helping identify editing sites that may not be detected through conventional de novo approaches.

5. **Post-recovery Integration and Visualization:** Integrates newly recovered sites into the orthology-based results and visualizes the updated conserved-site groups.

6. **Differential RNA Editing Analysis:** Identifies differences in RNA editing levels between experimental conditions or sample groups, supporting the investigation of condition-associated RNA editing patterns.

7. **Genomic Base Extraction:** Extracts nucleotide counts (A, C, G, and T) at specified genomic positions across samples, providing base-level information for downstream analyses and validation.

8. **Utility Scripts:** Provides additional tools for data preparation, coordinate conversion, site annotation, result integration, and other auxiliary operations required during RNA editing analyses.

## 3. Installation and Setup

### 3.1. Download KOBE

Download the KOBE software package and its Apptainer container image.

- **KOBE software:** `KOBE_v1.0.0.zip`
- **Apptainer container image:** Download link will be provided when the image is publicly available.

The Quick Test dataset is provided separately as `KOBE_Quick_Test_Data_v1.0.0.zip` and is also available in the repository under `test_data/`.

[GitHub Releases](https://github.com/mao1qaz/KOBE/releases)

### 3.2. Install KOBE

Extract the downloaded software package and run the installation script.

```bash
mkdir -p /path/to/software
cd /path/to/software

# Extract the software package
mkdir -p KOBE_v1.0.0
unzip -o KOBE_v1.0.0.zip -d KOBE_v1.0.0

# Install KOBE
bash KOBE_v1.0.0/KOBE/install.sh
```

Replace `/path/to/software` with your preferred installation directory.

### 3.3. Install Apptainer

KOBE uses Apptainer to run its containerized analysis environment. Apptainer 1.5.4 is recommended.

Create and activate a dedicated Conda environment:

```bash
# Create the Conda environment
conda create -y \
    -n apptainer-1.5.4 \
    --override-channels \
    -c conda-forge \
    apptainer=1.5.4 \
    python=3.11

# Activate the environment
conda activate apptainer-1.5.4
```

The Apptainer environment should be activated before running KOBE commands that require the container.

Download the KOBE container image and keep it in an accessible directory.

### 3.4. PATH Configuration

The installation script automatically creates symbolic links for all KOBE commands in `$HOME/.local/bin/` by default.

After creating the command links, the installer checks whether the installation directory is included in the system `PATH`.

- If the directory is already in `PATH`, no additional configuration is required.
- If the directory is not in `PATH`, the installer displays the required configuration command.

In the latter case, add the displayed command to `~/.bashrc` and reload the shell configuration:

```bash
source ~/.bashrc
```

After configuration, all KOBE commands can be executed directly from any working directory.

**Note:** The default command installation directory can be customized using the `RNA_EDITING_BIN_DIR` environment variable.

### 3.5. Verify Installation

**1. Check Apptainer installation**

```bash
apptainer --version
```

Expected output:

```text
apptainer version 1.5.4
```

**2. Check KOBE command availability**

```bash
command -v KOBE
```

Expected output (default installation):

```text
/home/username/.local/bin/KOBE
```

The displayed path depends on the command installation directory selected during installation.

**3. Check KOBE version**

```bash
KOBE --version
```

Expected output:

```text
KOBE 1.0.0
```

**4. Check container availability and write permissions**

Specify the paths to the downloaded container image and the project output directory:

```bash
SIF="/path/to/rna_editing.sif"
TEST_DIR="/path/to/project/output"
```

Run the following test:

```bash
mkdir -p "$TEST_DIR" &&
apptainer exec \
    --bind "$TEST_DIR:/tmp" \
    "$SIF" \
    sh -c 'f=$(mktemp /tmp/.kobe_write_test.XXXXXX) && rm "$f"' &&
echo "PASS: Container is accessible and the output directory is writable."
```

Expected output:

```text
PASS: Container is accessible and the output directory is writable.
```

Replace the example paths with the actual container image and project output directory. The test creates and removes a temporary file without modifying existing project data.

If all four checks pass, KOBE and its container environment are ready to use.

### 3.6. Running Long-running Commands in the Background

For long-running analyses, KOBE commands can be executed in the background to prevent interruption when the terminal session is closed.

**Recommended command format:**

```bash
nohup setsid --fork --wait KOBE <command> \
    [options] \
    > "/path/to/command.log" 2>&1 < /dev/null &
```

**Example:**

```bash
nohup setsid --fork --wait KOBE denovo_detect_map \
    --samples "/path/to/denovo_detect_data/samples.tsv" \
    --ref-dir "/path/to/denovo_detect_data/ref" \
    --out-dir "/path/to/denovo_detect_result" \
    --cores 30 \
    --sif "/path/to/rna_editing.sif" \
    > "/path/to/denovo_detect_map.log" 2>&1 < /dev/null &
```

The log file records standard output and error messages.

## 4. Quick Start

This example demonstrates orthology-based RNA editing analysis and conservation visualization using the two-species test dataset provided in `test_data/`.

**Prerequisites:** KOBE is installed, Apptainer is activated, and the KOBE container image is available.

Run the following commands from the project directory containing `test_data/`.

**1. Set paths**

```bash
PROJECT_DIR="$(pwd)"
TEST_DATA="$PROJECT_DIR/test_data"
OUT_DIR="$PROJECT_DIR/quick_test_output"
SIF="/path/to/rna_editing_fonts.sif"

mkdir -p "$OUT_DIR"
```

Replace the example SIF path with the actual container image location.

**2. Run orthology analysis**

```bash
KOBE ortho_edit \
    --cds_dir "$TEST_DATA/ortho_edit_data" \
    --edit_txt_dir "$TEST_DATA/ortho_edit_txt" \
    --work_dir "$OUT_DIR/ortho_result" \
    --threads 20 \
    --sif "$SIF" \
    --bind "$PROJECT_DIR"
```

**3. Visualize cross-species conservation**

```bash
KOBE conserve_plot \
    --annotated "$OUT_DIR/ortho_result/ortho_edit_result/merged_edit_sites_codon_table_annotated.csv" \
    --cds_dir "$TEST_DATA/ortho_edit_data" \
    --min_species 2 \
    --species_order "melanogaster,echinatior" \
    --source_priority "melanogaster,echinatior" \
    --outdir "$OUT_DIR/conserve_plot_result" \
    --sif "$SIF" \
    --col 8 \
    --bind "$PROJECT_DIR"
```

The results are saved in `quick_test_output/`.

For detailed instructions, parameter descriptions, and input/output formats, see the [KOBE User Manual](docs/KOBE_User_Manual.md).

## 5. Citation

If you use KOBE in your research, please cite:

*Citation information will be provided upon publication.*

## 6. License

KOBE is released under the [MIT License](LICENSE).
