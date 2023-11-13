# -*- coding: utf-8 -*-

import requests
import pandas as pd

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning) #Note genes are retrieved from mistDB's API using requests.get with disabled SSL verification, the mib-jouline-db.asc.ohio-state.edu server's certificate chain is incomplete according to SSL Labs. 

import os  # required for making the directory into which to save the plots of the TCS.

def retrieve_genes(genome_accession, page_max, sig_types):

    # get all genes of genome

    GENES = []
    counter = 0
    for pagenumber in range(1, page_max):
        url_genes = "https://mib-jouline-db.asc.ohio-state.edu/v1/genomes/{accession}/genes?page={page}&per_page=100".format(accession=genome_accession, page=pagenumber)
        r = requests.get(url_genes, verify=False) # Disabled SSL verification 
        json_data = r.json()
        genes = [dict for dict in json_data]
        if len(genes) != 0:
            GENES.extend(genes)
        counter = counter + len(json_data)
        if counter % 1000 == 0:
            print("still working", counter, "genes extracted")
        if len(json_data) == 0:  # loop stops if there are no more records on the page
            print()
            break

    # get signalling genes

    SIG_GENES = []
    for pagenumber in range( 1, page_max ):  # looping required by MIST3 database. MIST REST API returns a max of 100 results per page
        url_signal = "https://mib-jouline-db.asc.ohio-state.edu/v1/genomes/{accession}/signal-genes?page={page}&per_page=100".format( accession=genome_accession, page=pagenumber )
        json_data = requests.get(url_signal, verify=False).json()
        #json_data = r.json()
        if sig_types == "OCS, TCS":
            sig_genes = [ dict for dict in json_data if dict["ranks"][0] == "tcp" or dict["ranks"][0] == "ocp" ]  # list of both TCS and OCS genes
        elif sig_types == "OCS":
            sig_genes = [ dict for dict in json_data if dict["ranks"][0] == "ocp" ]  # list of just the OCS genes
        elif sig_types == "TCS":
            sig_genes = [ dict for dict in json_data if dict["ranks"][0] == "tcp" ]  # list of just the TCS genes
        if len(sig_genes) != 0:
            SIG_GENES.extend(sig_genes)
        if len(json_data) == 0:  # loop stops if there are no more records on the page
            break

    # find the additional info in the GENES list for any gene in the SIG_GENES list (using gene_id/id comparison) and then convert to a pandas dataframe
    SIG_GENES_PD = [
        {**sig_gene, **entry}
        for entry in GENES
        for sig_gene in SIG_GENES
        if sig_gene["gene_id"] == entry["id"]]
    SIG_GENES_PD = pd.DataFrame(SIG_GENES_PD)
    SIG_GENES_PD = SIG_GENES_PD.sort_values( by=["start"] )  # sort signalling genes sequentially by their location in the genome.

    GENES_PD = pd.DataFrame(GENES)
    GENES_PD = GENES_PD.sort_values( by=["start"] )  # sort so all genes are sorted in the df sequentially by where they appear in the genome.
    GENES_PD.loc[GENES_PD.start > GENES_PD.stop, "stop"] = ( GENES_PD.start + GENES_PD.stop )  # genomes are circular -> so some end positions on the last gene may be smaller than the start positions, for simplicity & plotting their positions are added together (so stop position is technically outside the genome length)

    cond = GENES_PD["id"].isin(SIG_GENES_PD["id"])
    GENES_PD.drop( GENES_PD[cond].index, inplace=True )  # remove the signalling genes from the ALL_GENES dataframe.

    HK_RR_SIG_GENES_PD=SIG_GENES_PD[SIG_GENES_PD.apply(lambda x: x['ranks']==['tcp','hk']or x['ranks'] == ['tcp', 'rr'], axis=1)]

    print( len(SIG_GENES_PD), sig_types, "HK, RR, HHK, HRR and other signalling genes were found in genome", genome_accession)
    print( len(HK_RR_SIG_GENES_PD), sig_types, " HK and RR only signalling genes were found in genome", genome_accession)


    return (HK_RR_SIG_GENES_PD, SIG_GENES_PD, GENES_PD)  # returns two data frames, SIG_GENES_PD is signalling genes, GENES_PD is all other genes.

def pair_TCS(SIG_GENES):

    # takes dataframe of all signalling genes and tries to find TCS pairs (HK & RR)

    TCS_GENES = SIG_GENES.loc[SIG_GENES["ranks"].apply(lambda x: len(x)) != 1]  # subset the signalling genes to only include TCS

    i = 0
    close_pairs = pd.DataFrame()
    correct_pairs = pd.DataFrame()
    strand_correct_pairs = pd.DataFrame()
    pairs_list = []
    orphans = pd.DataFrame()

    while i + 1 < len(TCS_GENES):
        gene1 = TCS_GENES.iloc[[i]]
        gene2 = TCS_GENES.iloc[[i + 1]]
        i += 1

        # Proximity (closer than 200bp)
        if abs(gene1.iloc[0]["stop"] - gene2.iloc[0]["start"]) <= 200:
            close_pairs = pd.concat([close_pairs,gene1,gene2],ignore_index=True)
            i+=1

            # Pair Type (HK & RR)
            if (gene1.iloc[0]["ranks"][1] == "hk" and gene2.iloc[0]["ranks"][1] == "rr") or (gene1.iloc[0]["ranks"][1] == "rr" and gene2.iloc[0]["ranks"][1] == "hk"):
                correct_pairs = pd.concat([correct_pairs,gene1, gene2],ignore_index=True)

                # Orientation (HK & RR are on same strand)
                if gene1.iloc[0]["strand"] == gene2.iloc[0]["strand"]:
                    strand_correct_pairs = pd.concat([strand_correct_pairs,gene1,gene2],ignore_index=True)
                    pairs_list.append(pd.concat([gene1, gene2],ignore_index=True))

    # Else TCS gene is an orphan
        else:
            orphans = pd.concat([orphans,gene1])

#making sure orphans and pairs match to total number of TCS genes, previous while loop misses a last orphan, so if TCS genes number doesn't match, append the last TCS gene which is an orphan.
        if len(strand_correct_pairs)+len(orphans) != len(TCS_GENES):
            orphans = pd.concat([orphans,TCS_GENES.iloc[[-1]]])


    print(
        "I found",
        int(len(close_pairs) / 2),
        "proximity pairs,",
        int(len(correct_pairs) / 2),
        "proximity & type pairs, ",
        int(len(strand_correct_pairs) / 2),
        "proximity & type & orientation correct pairs and",
        len(orphans),
        "orphans",
    )

    return pairs_list  # returns TCS pairs as a list, each element is dataframe containing information on the HK & partner RR

def get_TCS_target(TCS, ALL_GENES):

    # use extracted genes (ALL_GENES) to find what gene is closest to the TCS (that might be the target)
    
    canonical_orientation = 0
    
    canonical_TCS = []
    
    complete_TCS = []
    for pair in TCS:
        start_loc = pair.iloc[0]["start"]  # get the start location of the HK&RR pair
        end_loc = pair.iloc[1]["stop"]  # get the stop location of the HK&RR pair
        direct_loc = pair.iloc[0]["strand"]  # get the direction of the HK&RR pair
    
        left_neighbour = ALL_GENES[ALL_GENES.stop < start_loc].iloc[[-1]]  # get left and right neighbours of the HK&RR pair
        right_neighbour = ALL_GENES[ALL_GENES.start > end_loc].iloc[[0]]
    
        if (direct_loc == "+"):  # if the direction is + get left neighbour, else get right neighbouring gene.
            new = pd.concat([pair,left_neighbour],ignore_index=True)  # add the target gene to the dataframe containing the HK&RR
            complete_TCS.append(new)  # add the HK&RR&target df to the array of TCS
    
            if left_neighbour.iloc[0]["strand"] != "+":
                canonical_orientation += 1
                canonical_TCS.append(new)
    
        else:
            new = pd.concat([pair,right_neighbour],ignore_index=True)  # add the target gene to the dataframe containing the HK&RR
            complete_TCS.append(new)  # add the HK&RR&target df to the array of TCS
            
            if right_neighbour.iloc[0]["strand"] == "+":
                canonical_orientation += 1
                canonical_TCS.append(new)
    
    print(
        canonical_orientation,
        "TCS target pairs are in the canonical orientation, out of ",
        len(TCS),
        "total TCS systems",
    )

    return complete_TCS, canonical_TCS

def plot_TCS_systems(PAIRS_TARGETS_TCS, genome_accession, nickname):

    if not os.path.exists(str(nickname + " " + genome_accession)):
        os.makedirs(
            str(nickname + " " + genome_accession)
        )  # make directory into which to store gene plots (unless it already exists)
        os.makedirs(str(nickname + " " + genome_accession) + "/TCS_plots")

    for k in range(
        len(PAIRS_TARGETS_TCS)
    ):  # loop to get start, stop, strand, colours and names of each gene (HK, RR and target)
        starts = []
        stops = []
        strands = []
        colours = []
        names = []

        for i in range(3):

            # getting gene prodcut & name and using it as label for plot

            gene_product = PAIRS_TARGETS_TCS[k]["product"].iloc[i]
            if (gene_product is not None):  # removing superfluous info included in product description (to shorten label, ex. we know RR are proteins... etc.)
                for r in (
                    ("DNA-binding", ""),
                    ("sensor", ""),
                    ("protein", ""),
                    ("system", ""),
                    ("response regulator", "RR"),
                    ("two-component", ""),
                    ("histidine kinase", "HK"),
                ):
                    gene_product = gene_product.replace(*r)

            gene_name = str(PAIRS_TARGETS_TCS[k]["names"].iloc[i])  # get any additional name in 'names' column of the df
            final_name = (gene_product if gene_product is not None else "") + (gene_name if gene_name != "None" else "")
            names.append(final_name)

        for i in range(len(PAIRS_TARGETS_TCS[k])):
            start = PAIRS_TARGETS_TCS[k]["start"].iloc[i]
            starts.append(int(start))
            stop = PAIRS_TARGETS_TCS[k]["stop"].iloc[i]
            stops.append(int(stop))

            strand = PAIRS_TARGETS_TCS[k]["strand"].iloc[i]
            strands.append(strand)
            if i <= 1:
                type_TCS_gene = PAIRS_TARGETS_TCS[k]["ranks"].iloc[i][1]
                if type_TCS_gene == "hk":
                    colours.append("red")  # HK are red
                else:
                    colours.append("purple")  # RR are purple

        colours.append("green")  # target genes are green

        for n, i in enumerate(
            strands
        ):  # changing strand lists so it's in a format compatible with gds_features.add_feature
            if i == "-":
                strands[n] = -1
            elif i == "+":
                strands[n] = +1

        #####plotting####
        gdd = GenomeDiagram.Diagram("Test Diagram")
        gdt_features = gdd.new_track(
            1,  # make track onto which features (genes) will be plotted
            height=0.9,
            greytrack=True,
            greytrack_labels=0,
            # scale_smalltick_interval = 500,
            scale_largetick_interval=1000,
            scale_largeticks=0.25,
            scale_largetick_labels=0,
        )
        gds_features = gdt_features.new_set()

        for i in range(3):
            feature = SeqFeature(
                FeatureLocation(starts[i], stops[i]), strand=strands[i]
            )
            gds_features.add_feature(
                feature,
                sigil="ARROW",
                name=names[i],
                label=True,
                label_size=10,
                label_angle=15,
                colour=colours[i],
                arrowshaft_height=0.8,
            )

        gdd.draw(
            format="linear",
            pagesize=(25 * cm, 10 * cm),
            x=0.2,
            y=0.2,
            track_size=0.5,
            fragments=1,
            start=min(starts) - 100,
            end=max(stops) + 100,
        )

        file_name = (
            nickname
            + " "
            + genome_accession
            + "/TCS_plots"
            + "/TCS_system_"
            + str(k)
            + ".pdf"
        )
        gdd.write(file_name, "PDF")
    return

def get_sequences(PAIRS_TARGETS_TCS, genome_accession, left_spacer, right_spacer, nickname, sig_types):

    if not os.path.exists(str(nickname + " " + genome_accession)):
        os.makedirs(
            str(nickname + " " + genome_accession)
        )  # make directory into which to store gene plots (unless it already exists)
        os.makedirs(str(nickname + " " + genome_accession) + "/TCS_plots")


    # change the genome accesion (ga) number into a format that's compatible with ncbi's ftp folder organisation in order to find ftp directory (they split the number into groups of 3 numbers)
    ga = genome_accession.replace("_", "")
    split_ga = [ga[i : i + 3] for i in range(0, len(ga), 3)]
    ftp_target = ""
    for i in range(len(split_ga) - 1):
        ftp_target = ftp_target + split_ga[i] + "/"

    # download the gene annotations and dna from ncbi's ftp server for the specific ga.
    import ftputil  # from https://ftputil.sschwarzer.net/trac/wiki/Documentation#ftphost-objects

    host = ftputil.FTPHost(
        "ftp.ncbi.nlm.nih.gov", "anonymous", "password"
    )  # sign in to ncbi's ftp server.
    host.chdir(
        "/genomes/all/" + ftp_target
    )  # go to the genom accession number in the directory of the ftp
    ga_version = host.listdir(".")[
        -1
    ]  # get the latest version of the genome accession number files (if there are more than one)
    host.chdir(
        "/genomes/all/" + ftp_target + ga_version
    )  # go to the directory with the most recent version of the files.

    target_gff_file = (
        nickname + " " + genome_accession + "/" + genome_accession + ".gff.gz"
    )
    host.download(
        source=ga_version + "_genomic.gff.gz", target=target_gff_file  # downloads gff
    )  # saves gff as a target file

    target_fna_file = (
        nickname + " " + genome_accession + "/" + genome_accession + ".fna.gz"
    )
    host.download(
        source=ga_version + "_genomic.fna.gz", target=target_fna_file  # downloads fna
    )  # saves fna as a target file

    # unzips the files
    import gzip
    import shutil

    with gzip.open(target_fna_file, "rb") as f_in:
        with open(target_fna_file.replace("gz", ""), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    with gzip.open(target_gff_file, "rb") as f_in:
        with open(target_gff_file.replace("gz", ""), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Parse the annotations and sequence into Biopython format
    from BCBio import GFF
    from Bio import SeqIO

    # use Biopython to parse fna file
    in_seq_file = target_fna_file.replace("gz", "")
    in_seq_handle = open(in_seq_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    target_rec = []
    in_file = target_gff_file.replace("gz", "")
    in_handle = open(in_file)
    for rec in GFF.parse(in_handle, target_lines=30):  # adds whole fna to gff dictionary of features.
        target_rec.append(rec)
    in_handle.close()

    if len(target_rec) > 1:
        print(
            "Warning: there might be plasmid(s) in the gff file, checking that the sequence-region ID matches"
        )

    # get the record with the most features (genomic DNA)
    rec_features_len = []
    for rec in target_rec:
        rec_features_len.append(len(rec.features))
    rec = target_rec[rec_features_len.index(max(rec_features_len))]
    rec.seq = seq_dict[rec.id]

    # match the protein ids from MISTdb to the ones from ncbi ftp extracted sequences
    if nickname == 'EHEC':

        no_match_counter = 0
        for a in range(len(PAIRS_TARGETS_TCS)):
            TCS_group_df = PAIRS_TARGETS_TCS[a]
            fts = []

            offset = []

            for b in range(len(TCS_group_df)):

                mistdb_start = TCS_group_df.start.iloc[b]
                mistdb_end = TCS_group_df.stop.iloc[b]
                mistdb_strand = TCS_group_df.strand.iloc[b]

                ft_seq = []

                RefSeq_ft_id = "gene" + TCS_group_df.stable_id.iloc[b].replace(
                    genome_accession, ""
                )
                for i in range(len(rec.features)):
                    ft = rec.features[
                        i
                    ]  # search through the NCBI downloaded genome features to find matching feature
                    if ft.id.replace('_','') == RefSeq_ft_id:  # if you find matching feature id
                        # append the sequence, from start and end sites directly in ncbi
                        ft_seq = rec.seq[
                            (ft.location.start - left_spacer) : (
                                ft.location.end + right_spacer
                            )
                        ].seq
                        if ft.strand == 1:
                            ft_seq = ft_seq
                        else:
                            ft_seq = (
                                ft_seq.reverse_complement()
                            )  # get reverse complement of the sequence for features on the reverse strand
                        fts.append(str(ft_seq))

                        # if the start and end sites for matched ids are not the same, record this.
                        if ft.location.start == mistdb_start:
                            print("mistdb and ncbi are cohesive")
                        else:
                            offset.append(ft.location.start - mistdb_start)

                if (
                    ft_seq == []
                ):  # if we don't find a feature in ncbi features with matching stable id, use the coordinates from MISTdb

                    no_match_counter = no_match_counter + 1

                    empty = []
                    # double check if there is an offset between MISTdb coordinates and ncbi feature coordinates, if there is, take this into account.
                    if offset != empty:
                        ft_seq = rec.seq[
                            (mistdb_start + offset[0] - left_spacer) : (
                                mistdb_end + offset[0] + right_spacer
                            )
                        ].seq
                    else:
                        ft_seq = rec.seq[
                            (mistdb_start - left_spacer) : (mistdb_end + right_spacer)
                        ].seq

                    # append the correct sequence (depending on strand)
                    if mistdb_strand == "+":
                        ft_seq = ft_seq
                    else:
                        ft_seq = (
                            ft_seq.reverse_complement()
                        )  # get reverse complement of the sequence for features on the reverse strand
                    fts.append(str(ft_seq))

            TCS_group_df["sequence"] = fts
        
    else:
        no_match_counter = 0
        for a in range(len(PAIRS_TARGETS_TCS)):
            TCS_group_df = PAIRS_TARGETS_TCS[a]
            fts = []

            offset = []

            for b in range(len(TCS_group_df)):

                mistdb_start = TCS_group_df.start.iloc[b]
                mistdb_end = TCS_group_df.stop.iloc[b]
                mistdb_strand = TCS_group_df.strand.iloc[b]

                ft_seq = []

                RefSeq_ft_id = "gene" + TCS_group_df.stable_id.iloc[b].replace(
                    genome_accession, ""
                )
                for i in range(len(rec.features)):
                    ft = rec.features[
                        i
                    ]  # search through the NCBI downloaded genome features to find matching feature
                    if ft.id == RefSeq_ft_id:  # if you find matching feature id
                        # append the sequence, from start and end sites directly in ncbi
                        ft_seq = rec.seq[
                            (ft.location.start - left_spacer) : (
                                ft.location.end + right_spacer
                            )
                        ].seq
                        if ft.strand == 1:
                            ft_seq = ft_seq
                        else:
                            ft_seq = (
                                ft_seq.reverse_complement()
                            )  # get reverse complement of the sequence for features on the reverse strand
                        fts.append(str(ft_seq))

                        # if the start and end sites for matched ids are not the same, record this.
                        if ft.location.start == mistdb_start:
                            print("mistdb and ncbi are cohesive")
                        else:
                            offset.append(ft.location.start - mistdb_start)

                if (
                    ft_seq == []
                ):  # if we don't find a feature in ncbi features with matching stable id, use the coordinates from MISTdb

                    no_match_counter = no_match_counter + 1

                    empty = []
                    # double check if there is an offset between MISTdb coordinates and ncbi feature coordinates, if there is, take this into account.
                    if offset != empty:
                        ft_seq = rec.seq[
                            (mistdb_start + offset[0] - left_spacer) : (
                                mistdb_end + offset[0] + right_spacer
                            )
                        ].seq
                    else:
                        ft_seq = rec.seq[
                            (mistdb_start - left_spacer) : (mistdb_end + right_spacer)
                        ].seq

                    # append the correct sequence (depending on strand)
                    if mistdb_strand == "+":
                        ft_seq = ft_seq
                    else:
                        ft_seq = (
                            ft_seq.reverse_complement()
                        )  # get reverse complement of the sequence for features on the reverse strand
                    fts.append(str(ft_seq))

            TCS_group_df["sequence"] = fts

    pd.concat(PAIRS_TARGETS_TCS, keys=list(range(len(PAIRS_TARGETS_TCS)))).to_csv(
        nickname + " " + genome_accession + "/paired_TCS_targets_and_Sequences.csv"
    )
    finish_line = "Your " + sig_types + " targets and sequences have been saved"
    print(finish_line)
    print(
        "For ",
        no_match_counter,
        " features I found no ID matches for Mistdb stable IDs in ncbi features, so used MistDB start and stop positions instead",
    )
    return

def get_combined_sequences(PAIRS_TARGETS_TCS, genome_accession, left_spacer, right_spacer, nickname, sig_types):

    # change the genome accesion (ga) number into a format that's compatible with ncbi's ftp folder organisation in order to find ftp directory (they split the number into groups of 3 numbers)
    ga = genome_accession.replace("_", "")
    split_ga = [ga[i : i + 3] for i in range(0, len(ga), 3)]
    ftp_target = ""
    for i in range(len(split_ga) - 1):
        ftp_target = ftp_target + split_ga[i] + "/"

    # download the gene annotations and dna from ncbi's ftp server for the specific ga.
    import ftputil  # from https://ftputil.sschwarzer.net/trac/wiki/Documentation#ftphost-objects

    host = ftputil.FTPHost("ftp.ncbi.nlm.nih.gov", "anonymous", "password")  # sign in to ncbi's ftp server.
    host.chdir("/genomes/all/" + ftp_target)  # go to the genom accession number in the directory of the ftp
    ga_version = host.listdir(".")[-1]  # get the latest version of the genome accession number files (if there are more than one)
    host.chdir("/genomes/all/" + ftp_target + ga_version)  # go to the directory with the most recent version of the files.

    target_gff_file = (nickname + " " + genome_accession + "/" + genome_accession + ".gff.gz")
    host.download(source=ga_version + "_genomic.gff.gz", target=target_gff_file)  # saves gff as a target file # downloads gff

    target_fna_file = (nickname + " " + genome_accession + "/" + genome_accession + ".fna.gz")
    host.download(source=ga_version + "_genomic.fna.gz", target=target_fna_file)  # saves fna as a target file # downloads fna

    # unzips the fna (DNA) and gff (annotation) files
    import gzip
    import shutil

    with gzip.open(target_fna_file, "rb") as f_in:
        with open(target_fna_file.replace("gz", ""), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    with gzip.open(target_gff_file, "rb") as f_in:
        with open(target_gff_file.replace("gz", ""), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    # Parse the annotations and sequence into Biopython format
    from BCBio import GFF
    from Bio import SeqIO

    # use Biopython to parse fna file
    in_seq_file = target_fna_file.replace("gz", "")
    in_seq_handle = open(in_seq_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    target_rec = []
    in_file = target_gff_file.replace("gz", "")
    in_handle = open(in_file)
    for rec in GFF.parse(in_handle, target_lines=30):  # adds whole fna to gff dictionary of features.
        target_rec.append(rec)
    in_handle.close()

    if len(target_rec) > 1:
        print("Warning: there might be plasmid(s) in the gff file, checking that the sequence-region ID matches")

    # get the record with the most features (probably the genomic)
    rec_features_len = []
    for rec in target_rec:
        rec_features_len.append(len(rec.features))
    rec = target_rec[rec_features_len.index(max(rec_features_len))]
    rec.seq = seq_dict[rec.id]

    # match the protein ids from MISTdb to the ones from ncbi ftp extracted sequences
    if nickname == 'EHEC':

        no_match_counter = 0
        for a in range(len(PAIRS_TARGETS_TCS)):
            TCS_group_df = PAIRS_TARGETS_TCS[a]
            fts = []

            offset = []

            for b in range(len(TCS_group_df)):

                mistdb_start = TCS_group_df.start.iloc[b]
                mistdb_end = TCS_group_df.stop.iloc[b]
                mistdb_strand = TCS_group_df.strand.iloc[b]

                ft_seq = []

                RefSeq_ft_id = "gene" + TCS_group_df.stable_id.iloc[b].replace(
                    genome_accession, ""
                )
                for i in range(len(rec.features)):
                    ft = rec.features[
                        i
                    ]  # search through the NCBI downloaded genome features to find matching feature
                    if ft.id.replace('_','') == RefSeq_ft_id:  # if you find matching feature id
                        # append the sequence, from start and end sites directly in ncbi
                        ft_seq = rec.seq[
                            (ft.location.start - left_spacer) : (
                                ft.location.end + right_spacer
                            )
                        ].seq
                        if ft.strand == 1:
                            ft_seq = ft_seq
                        else:
                            ft_seq = (
                                ft_seq.reverse_complement()
                            )  # get reverse complement of the sequence for features on the reverse strand
                        fts.append(str(ft_seq))

                        # if the start and end sites for matched ids are not the same, record this.
                        if ft.location.start == mistdb_start:
                            print("mistdb and ncbi are cohesive")
                        else:
                            offset.append(ft.location.start - mistdb_start)

                if (
                    ft_seq == []
                ):  # if we don't find a feature in ncbi features with matching stable id, use the coordinates from MISTdb

                    no_match_counter = no_match_counter + 1

                    empty = []
                    # double check if there is an offset between MISTdb coordinates and ncbi feature coordinates, if there is, take this into account.
                    if offset != empty:
                        ft_seq = rec.seq[
                            (mistdb_start + offset[0] - left_spacer) : (
                                mistdb_end + offset[0] + right_spacer
                            )
                        ].seq
                    else:
                        ft_seq = rec.seq[
                            (mistdb_start - left_spacer) : (mistdb_end + right_spacer)
                        ].seq

                    # append the correct sequence (depending on strand)
                    if mistdb_strand == "+":
                        ft_seq = ft_seq
                    else:
                        ft_seq = (
                            ft_seq.reverse_complement()
                        )  # get reverse complement of the sequence for features on the reverse strand
                    fts.append(str(ft_seq))

            TCS_group_df["sequence"] = fts
        
    else:
        no_match_counter = 0

        canonical_fts = pd.DataFrame(data={'Name':[], 'Sequence':[]}) #list of combined TCS pair sequences. 
        group_ft_seqs = []

        for a in range(len(PAIRS_TARGETS_TCS)): #for each group of TCS
            TCS_group_df = PAIRS_TARGETS_TCS[a]
            
            fts = [] #list of DNA sequences
            
            canonical_fts = pd.DataFrame(data={'Name':[], 'Sequence':[]}) #list of combined TCS pair sequences. 
            
            offset = []
            
            TCS1_start = TCS_group_df.start.iloc[0]
            TCS1_stop = TCS_group_df.stop.iloc[0]
            TCS2_start = TCS_group_df.start.iloc[1]
            TCS2_stop = TCS_group_df.stop.iloc[1]

            Target_start=TCS_group_df.start.iloc[2]
            Target_stop=TCS_group_df.stop.iloc[2]

            if Target_stop < TCS1_start:
                TCS1_start = Target_stop
            else: 
                TCS2_stop = Target_start
            
            group_seq = rec.seq[(TCS1_start - left_spacer) : (TCS2_stop + right_spacer)].seq
            group_ft_seqs.append(str(group_seq))
            

            for b in range(len(TCS_group_df)):

                mistdb_start = TCS_group_df.start.iloc[b]
                mistdb_end = TCS_group_df.stop.iloc[b]
                mistdb_strand = TCS_group_df.strand.iloc[b]

                ft_seq = []

                RefSeq_ft_id = "gene" + TCS_group_df.stable_id.iloc[b].replace(genome_accession, "")
                for i in range(len(rec.features)):
                    ft = rec.features[i]  # search through the NCBI downloaded genome features to find matching feature
                    
                    if ft.id == RefSeq_ft_id:  # if you find matching feature id
                        
                        # append the sequence, from start and end sites directly in ncbi
                        ft_seq = rec.seq[(ft.location.start - left_spacer) : (ft.location.end + right_spacer)].seq
                        if ft.strand == 1:
                            ft_seq = ft_seq
                        else:
                            ft_seq = (ft_seq.reverse_complement())  # get reverse complement of the sequence for features on the reverse strand
                        fts.append(str(ft_seq))
                        

                        # if the start and end sites for matched ids are not the same, record this.
                        if ft.location.start == mistdb_start:
                            print("mistdb and ncbi are cohesive")
                        else:
                            offset.append(ft.location.start - mistdb_start)

                if (ft_seq == []):  # if we don't find a feature in ncbi features with matching stable id, use the coordinates from MISTdb

                    no_match_counter = no_match_counter + 1

                    empty = []
                    # double check if there is an offset between MISTdb coordinates and ncbi feature coordinates, if there is, take this into account.
                    if offset != empty:
                        ft_seq = rec.seq[(mistdb_start + offset[0] - left_spacer) : (mistdb_end + offset[0] + right_spacer)].seq
                    else:
                        ft_seq = rec.seq[(mistdb_start - left_spacer) : (mistdb_end + right_spacer)].seq

                    # append the correct sequence (depending on strand)
                    if mistdb_strand == "+":
                        ft_seq = ft_seq
                    else:
                        ft_seq = (ft_seq.reverse_complement())  # get reverse complement of the sequence for features on the reverse strand
                    fts.append(str(ft_seq))
            TCS_group_df["sequence"] = fts
        canonical_fts["Sequence"] = group_ft_seqs
        canonical_fts["Name"] = get_simplified_gene_names(PAIRS_TARGETS_TCS)

    pd.concat(PAIRS_TARGETS_TCS, keys=list(range(len(PAIRS_TARGETS_TCS)))).to_csv(nickname + " " + genome_accession + "/paired_TCS_targets_and_Sequences.csv")
    
    #count SapI & BsaI sites & add columns with count of RE sites. 
    added_sites_seq_pd, SapI_compatible_TCS, BsaI_compatible_TCS, BsaI_not_SapI_compatible_TCS, RE_sites_df = count_BsaI_and_SapI_sites(canonical_fts)
    
    #save tables to csv
    BsaI_compatible_TCS.to_csv(nickname + " " + genome_accession + "/" + nickname + "BsaI_compatible_canonical_TCS.csv")
    SapI_compatible_TCS.to_csv(nickname + " " + genome_accession + "/" + nickname + "SapI_compatible_canonical_TCS.csv")
    added_sites_seq_pd.to_csv(nickname + " " + genome_accession + "/" + nickname + "canonical_grouped_seq.csv")
    BsaI_not_SapI_compatible_TCS.to_csv(nickname + " " + genome_accession + "/" + nickname + "BsaI_not_SapI_compatible_TCS.csv")
    RE_sites_df.to_csv(nickname + " " + genome_accession + "/" + nickname + "RE_sites.csv")

    finish_line = "Your " + sig_types + " targets and sequences have been saved"
    print(finish_line)
    print("For ",no_match_counter," features I found no ID matches for Mistdb stable IDs in ncbi features, so used MistDB start and stop positions instead")
    
    return canonical_fts, RE_sites_df

def complete_sensor_search(genome_accession, nickname, page_max, sig_types, left_spacer, right_spacer):
    # function that runs the complete sensor search pipeline & saves all the recorded information into a folder

    import time

    start = time.time()

    SIG_GENES, ALL_GENES = retrieve_genes(genome_accession, page_max, sig_types)  # this retrieves all genes (signalling and none) from the MIST database
    PAIRS_TCS = pair_TCS(SIG_GENES)  # this finds TCS pairs (HK&RR)
    PAIRS_TARGETS_TCS, CANONICAL_PAIRS_TARGETS = get_TCS_target(PAIRS_TCS, ALL_GENES)  # this find the TCS target gene
    #plot_TCS_systems(
    #    PAIRS_TARGETS_TCS, genome_accession, nickname
    #)  # this plots the HK&RR & target
    get_sequences(
        PAIRS_TARGETS_TCS,
        genome_accession,
        left_spacer,
        right_spacer,
        nickname,
        sig_types,
    )

    get_sequences(
        CANONICAL_PAIRS_TARGETS,
        genome_accession,
        left_spacer,
        right_spacer,
        nickname,
        sig_types,
    )

    end = time.time()
    print("I took", end - start, "seconds to run")

    return PAIRS_TARGETS_TCS, CANONICAL_PAIRS_TARGETS

def combined_seq_complete_sensor_search(genome_accession, nickname, page_max, sig_types, left_spacer, right_spacer):
    # function that runs the complete sensor search pipeline & saves all the recorded information into a folder

    import time

    start = time.time()

    SIG_GENES, ALL_GENES = retrieve_genes(genome_accession, page_max, sig_types)  # this retrieves all genes (signalling and none) from the MIST database
    PAIRS_TCS = pair_TCS(SIG_GENES)  # this finds TCS pairs (HK&RR)
    PAIRS_TARGETS_TCS, CANONICAL_PAIRS_TARGETS = get_TCS_target(PAIRS_TCS, ALL_GENES)  # this find the TCS target gene
    #plot_TCS_systems(
    #    PAIRS_TARGETS_TCS, genome_accession, nickname
    #)  # this plots the HK&RR & target
    get_combined_sequences(
        CANONICAL_PAIRS_TARGETS,
        genome_accession,
        left_spacer,
        right_spacer,
        nickname,
        sig_types,
    )

    end = time.time()
    print("I took", end - start, "seconds to run")

    return CANONICAL_PAIRS_TARGETS

def align(table, table2):

    from Bio.Seq import Seq
    from Bio import pairwise2
    import pandas as pd

    # adding amino acid sequences to tables
    for j in range(len(table)):
        protein_seqs = []
        for i in range(3):
            coding_dna = Seq(table[j].iloc[i]["sequence"])
            protein_seq = coding_dna.translate()
            protein_seqs.append(protein_seq)
        table[j]["aa_seq"] = protein_seqs
    for j in range(len(table2)):
        protein_seqs = []
        for i in range(3):
            coding_dna = Seq(table2[j].iloc[i]["sequence"])
            protein_seq = coding_dna.translate()
            protein_seqs.append(protein_seq)
        table2[j]["aa_seq"] = protein_seqs

    # aligning all  elements in TCS group of table1 to elements in TCS groupe of table2.

    all_alignments = []
    scores = []
    len_seq1 = []
    name_seq1 = []
    len_seq2 = []
    name_seq2 = []
    stableid_seq1 = []
    stableid_seq2 = []

    for j in range(len(table)):
        for m in range(2):
            seq1 = table[j].iloc[m]["aa_seq"]

            for k in range(len(table2)):
                for l in range(2):
                    seq2 = table2[k].iloc[l]["aa_seq"]

                    alignments = pairwise2.align.globalxx(seq1, seq2)
                    all_alignments.append(alignments)

                    group_scores = []
                    for alignment in alignments:
                        group_scores.append(alignment.score)
                    scores.append(max(group_scores))

                    len_seq1.append(len(seq1))
                    len_seq2.append(len(seq2))

                    name_seq1.append(
                        str(table[j].iloc[m]["product"])
                        + str(table[j].iloc[m]["names"])
                    )
                    name_seq2.append(
                        str(table2[k].iloc[l]["product"])
                        + str(table2[k].iloc[l]["names"])
                    )

                    stableid_seq1.append(table[j].iloc[m]["stable_id"])
                    stableid_seq2.append(table2[k].iloc[l]["stable_id"])

    dfscores = pd.DataFrame(
        data={
            "alignment_scores": scores,
            "seq1_len": len_seq1,
            "seq2_len": len_seq2,
            "name_seq1": name_seq1,
            "name_seq2": name_seq2,
            "stableid_seq1": stableid_seq1,
            "stableid_seq2": stableid_seq2,
        }
    )

    return dfscores

def plot_alignment_results(dfscores):
    sorted_by_score = dfscores.sort_values("alignment_scores", ascending=False)
    sorted_by_score["proportional"] = (
        sorted_by_score["alignment_scores"] / sorted_by_score["seq1_len"]
        + sorted_by_score["alignment_scores"] / sorted_by_score["seq2_len"]
    ) / 2
    sorted_by_proportional = sorted_by_score.sort_values(
        "proportional", ascending=False
    )
    sorted_by_proportional = sorted_by_proportional.reset_index().reset_index()
    align_plot = sorted_by_proportional.plot.scatter(
        x="level_0", y="proportional", c="black", s=1
    )
    return align_plot

def get_simplified_gene_names(PAIRS_TARGETS_TCS):
    names=[]
    for k in range(len(PAIRS_TARGETS_TCS)):
        group_name=''
        for i in range(2):
            gene_product = PAIRS_TARGETS_TCS[k]["product"].iloc[i]
            
            if (gene_product is not None):  # removing superfluous info included in product description (to shorten label, ex. we know RR are proteins... etc.)
                for r in (("DNA-binding", ""),
                            ("sensor", ""),
                            ("protein", ""),
                            ("system", ""),
                            ("response regulator", "RR"),
                            ("two-component", ""),
                            ("histidine kinase", "HK"),
                            ("transcription factor", "TF")):
                            gene_product = gene_product.replace(*r)

                gene_name = str(PAIRS_TARGETS_TCS[k]["names"].iloc[i])  # get any additional name in 'names' column of the df
                final_name = (gene_product if gene_product is not None else "") + (gene_name if gene_name != "None" else "")
                group_name = group_name+final_name

        names.append(group_name)
    return names

def count_BsaI_and_SapI_sites(seq_pd):
    BsaI_sites=[]
    SapI_sites=[]
    PaqCI_sites=[]
    BsaI_site='GGTCTC'
    BsaI_rev_site='GAGACC'
    SapI_site='GCTCTTC'
    SapI_rev_site='GAAGAGC'
    PaqCI_site='CACCTGC'
    PaqCI_rev_site='GTGGACG'
    
    for i in range(len(seq_pd)):
        comb_seq = seq_pd['Sequence'].iloc[i]
        N_BsaI_sites = comb_seq.count(BsaI_site)+comb_seq.count(BsaI_rev_site)
        N_SapI_sites = comb_seq.count(SapI_site)+comb_seq.count(SapI_rev_site)
        N_PaqCI_sites = comb_seq.count(PaqCI_site)+comb_seq.count(PaqCI_rev_site)
        BsaI_sites.append(N_BsaI_sites)
        SapI_sites.append(N_SapI_sites)
        PaqCI_sites.append(N_PaqCI_sites)
    seq_pd['BsaI_sites']=BsaI_sites
    seq_pd['SapI_sites']=SapI_sites
    seq_pd['PaqCI_sites']=PaqCI_sites

    SapI_compatible_TCS=seq_pd[(seq_pd['SapI_sites']==0)]
    BsaI_compatible_TCS=seq_pd[(seq_pd['BsaI_sites']==0)]
    BsaI_not_SapI_compatible_TCS=BsaI_compatible_TCS[(BsaI_compatible_TCS['SapI_sites']>0)]
    PaqCI_compatible_TCS=seq_pd[(seq_pd['PaqCI_sites']==0)]


    TCS_with_SapI_sites=seq_pd[(seq_pd['SapI_sites']>0)]
    TCS_with_BsaI_sites=seq_pd[(seq_pd['BsaI_sites']>0)]
    TCS_with_PaqCI_sites=seq_pd[(seq_pd['PaqCI_sites']>0)]

    TCS_with_BsaI_and_SapI_sites=seq_pd[(seq_pd['BsaI_sites']>0)&(seq_pd['SapI_sites']>0)]
    print('Out of a total', len(seq_pd), 'two canonical component systems')
    
    print(len(TCS_with_SapI_sites), 'have SapI sites')
    print(len(TCS_with_BsaI_sites), 'have BsaI sites')
    print(len(BsaI_not_SapI_compatible_TCS),'have SapI sites but no BsaI sites')
    print(len(TCS_with_PaqCI_sites), 'have PaqCI sites')


    d = {'Canonical_TCS': [len(seq_pd)],
     'SapI_compatible_TCS': [len(SapI_compatible_TCS)],
     'BsaI_compatible_TCS': [len(BsaI_compatible_TCS)],
     'BsaI_not_SapI_compatible_TCS': [len(BsaI_not_SapI_compatible_TCS)],
     'PaqCI_compatible_TCS':[len(PaqCI_compatible_TCS)]}
    RE_sites_df = pd.DataFrame(data=d)

    print(len(TCS_with_BsaI_and_SapI_sites), 'have both SapI and BsaI sites')
    print('Leaving',len(SapI_compatible_TCS), 'SapI compatible TCS, and',len(BsaI_compatible_TCS), 'BsaI compatible TCS, out of',len(seq_pd), 'total TCS')

    return seq_pd, SapI_compatible_TCS, BsaI_compatible_TCS, BsaI_not_SapI_compatible_TCS, RE_sites_df







