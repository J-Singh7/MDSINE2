# Replication of MDSINE2 results

### Internal documentation for manual steps (not necessary for running MDSINE2 or the tutorials)

0. Generate the `qpcr.tsv` table, `metadata.tsv` table, and `perturbations.tsv` table.
1. Run DADA2 to get ASVs. For people with access to ErisOne, this is in `/data/cctm/darpa_perturbation_mouse_study/dada_travis_out_rdp18`. Move the following files:
    * `/data/cctm/darpa_perturbation_mouse_study/dada_travis_out_rdp18/counts.tsv` --> `MDSINE2_Paper/datasets/gibson/counts.tsv`
    * `/data/cctm/darpa_perturbation_mouse_study/dada_travis_out_rdp18/rdp_species.tsv` --> `MDSINE2_Paper/datasets/gibson/rdp_species.tsv`
    * `/data/cctm/darpa_perturbation_mouse_study/dada_travis_out_rdp18/silva_species.tsv` --> `MDSINE2_Paper/datasets/gibson/silva_species.tsv`
    
2. Align and then place ASV sequences on reference phylogenetic tree. For people with access to ErisOne, this is in `/data/cctm/darpa_perturbation_mouse_study/phylogenetic_placement/run_phyloplacement_ASVs.sh`. Make sure you put the `fasta` file of the ASV sequences from DADA2 into the `phylogenetic_placement/data/query_reads/reads_ASV.fa`. Once finished running, move the file `phylogenetic_placement/output_ASVs/placed_seqs_on_v4_region.sto` into MDSINE2 folder `MDSINE2_Paper/gibson_dataset/files/phylogenetic_placement_ASVs/placed_seqs_on_v4_region.sto`. 
3. Agglomerate the ASVs into OTUs
    Once we have the sequences aligned, we can agglomerate the ASVs into OTUs. This is done with the script:
    ```bash
    ./preprocessing_agglomeration.sh
    ```
    This script agglomerates the ASVs into OTUs and generates consensus sequences. The above script writes to the folder `MDSINE2_Paper/processed_data`.
4. One of the outputs of `./preprocessing_agglomeration.sh` is a fasta file. Use any of the fasta files from `MDSINE2_Paper/processed_data` (`gibson_healthy_agg.fa`, `gibson_uc_agg.fa`, ...) for the input for phylogenetic placement. For people with access to ErisOne, move any file to `/data/cctm/darpa_perturbation_mouse_study/phylogenetic_placement/data/query_reads/reads_OTU.fa`. Once moved, run the script `/data/cctm/darpa_perturbation_mouse_study/phylogenetic_placement/run_phyloplacement_OTUs.sh`. This will align and place the OTU sequences on a phylogenetic tree. Move the following files:
    * `phylogenetic_placement/output_OTUs/newick_tree_full_taxid.nhx` --> `MDSINE2_Paper/gibson_dataset/files/phylogenetic_placement_OTUs/phylogenetic_tree_full_taxid.nhx`
        - This is a newick tree of the OTUs placed onto the reference tree. The leaf names are the OTU names and the **sequence ID** of the reference sequences
    * `phylogenetic_placement/output_OTUs/newick_tree_full_speciesName.nhx` --> `MDSINE2_Paper/gibson_dataset/files/phylogenetic_placement_OTUs/phylogenetic_tree_full.nhx`
        - This is a newick tree of the OTUs placed onto the reference tree. The leaf names are the OTU names and the **species names** of the reference sequences. This is the exact same tree as the prior one.
    * `phylogenetic_placement/output_OTUs/newick_tree_query_reads.nhx` --> `MDSINE2_Paper/gibson_dataset/files/phylogenetic_placement_OTUs/phylogenetic_tree_only_query.nhx`
        - This is a newick tree of only the OTUs of the previous tree. The only thing done are the reference sequences are pruned off the tree leaving only the OTUs. This is done for efficiency when making figures.
5. The last offline preprocessing is assigning taxonomy of the OTU consensus sequences with confidences with a naive Bayes classifier implemented with RDP. You can run this on the RDP website. Place that file in `MDSINE2_Paper/gibson_dataset/files/assign_taxonomy_OTUs/taxonomy_RDP.txt`. This file is used when generating the consensus taxonomy of the OTUs. This is used when the taxonomy of the ASVs do not agree with each other.