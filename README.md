Script outputs (e.g., HTML files, visualizations) are included here (contrary to standard practice) to support result sharing and transparency.

"To maintain code and outputs in an organised, accessible, and interpretable manner over time, a clear and consistent directory structure was adopted. All resources are stored under the root project directory (aef_forensic_metagenomics) as such:

 - data/ for static resources like adapters, indices, and databases (for trimming, alignment, and classification);
 - scripts/ for code, with helper_scripts/ for custom helper modules and logs/ for script outputs;
 - and tools/ for third-party software. 

Study-specific directories can be added alongside these, each containing a raw_data/ folder with read-only paired-end FASTQ files to preserve data integrity. Upon pipeline execution, processed_data/ and results/ folders are generated automatically. Each pipeline execution, referred to as “run”, receives a unique timestamp ID and is designated a subfolder within results/runs/. 

Files in processed_data/ are stored temporarily and overwritten by subsequent runs. 

This structure enables scripts to locate expected input and output files automatically."
