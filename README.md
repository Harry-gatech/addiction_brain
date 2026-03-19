Addiction Brain Analysis Pipeline
Overview
The addiction_brain package is a domain-specific research pipeline designed for neuroscience and chemoinformatics. It provides tools for analyzing drug-target interactions, pathway enrichment, and brain-region-specific data. The package is tailored for researchers working on addiction, brain disorders, and organ-level systems biology.

Features
Drug-Target Interaction Analysis: Tools for mapping and analyzing drug-target interactions using cheminformatics libraries like RDKit.
Pathway Enrichment Analysis: Automated workflows for identifying enriched pathways and performing statistical tests.
Brain-Region-Specific Analysis: Functions for analyzing targets and pathways specific to brain regions.
Data Deduplication: Utilities for cleaning and deduplicating chemical datasets.
Visualization: Support for generating figures and visual summaries of results.
Installation
To use the addiction_brain package, clone the repository and install the required dependencies:
# Clone the repository
git clone https://github.com/your-repo/addiction_brain.git

# Navigate to the project directory
cd addiction_brain

# Install dependencies
pip install -r requirements.txt

Usage
1. Deduplication
Use the deduplication module to clean and deduplicate chemical datasets:

from addiction_brain.deduplication import clean_inchi_butina

cleaned_data = clean_inchi_butina("path/to/mapping_file.txt")

2. BBB Permeability Analysis
Analyze blood-brain barrier (BBB) permeability using the BBB_calc module:

from addiction_brain.BBB_calc import annotate_bbb

annotated_data = annotate_bbb(dataframe, smiles_col="smiles")

3. Pathway Enrichment
Perform pathway enrichment analysis with the Enrichment_target module:

from addiction_brain.Enrichment_target import enrichment

enriched_data = enrichment(binding_data, N_add_total, N_non_total)

4. Network Analysis
Analyze drug trajectories and network interactions using the network_hop module:

from addiction_brain.network_hop import run_trajectory

results, summary = run_trajectory(network_df, add_targets, non_targets)


Project Structure
Dependencies
Python 3.8+
pandas
numpy
RDKit
scipy
statsmodels
tqdm


Contributing
Contributions are welcome! If you have ideas for new features or improvements, feel free to open an issue or submit a pull request.

License
This project is licensed under the MIT License. See the LICENSE file for details.

Acknowledgments
Special thanks to the contributors and the open-source community for providing the tools and libraries that made this project possible.

