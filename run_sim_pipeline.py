import sys
import os

for file in os.listdir(os.getcwd()):
    if os.path.isdir(file):
        sys.path.append(os.path.join(os.getcwd(), file))
src_path = os.getcwd()

# # the followings are simulation parameters

subclone_num = 1    # number of simulated subclones
subclone_ratio = "100"  # the corresponding subclone ratio (please refer to VISOR for its format)

dir = './'  # output dir
if not os.path.exists(dir):
    os.mkdir(dir)
    
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15' ,'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'] # simulated chroms


chr_sv_muns = {'chr1':800, 'chr2': 800, 'chr3': 700, 'chr4': 700, 'chr5': 600, 'chr6': 500, 'chr7': 500, 'chr8': 500, 'chr9': 500, 'chr10': 500, 'chr11': 500, 'chr12': 500, 'chr13': 400, 'chr14': 320, 'chr15': 320, 'chr16': 320, 'chr17': 300, 'chr18': 300, 'chr19': 200, 'chr20': 200, 'chr21': 200, 'chr22': 200, 'chrX': 500}   # number of simulated SV for each chrom

chr_snp_path = '../templatewithsnps/'   # template fasta file with inserted SNP (refer to VISOR for more)
ref_path = '/home/DATA/REFGENOMEDB/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa' # path of ref file
N_region_path = '../N_regions/' # path of N regions

# N_region_path = '/mnt/d/Workspace/test/exclude.bed'
read_type = 'PB'
if read_type not in ["ONT", "PB"]:
    print("Not a correct read type, must be ONT or PB!!!!!!!!!!")
    exit(0)

coverage = 5
threads = 3
sv_len_mean = 750   # SV mean length
sv_len_s = 150  # SV length deviation 
basic_sv_type_from_visor = 'insertion,deletion,inversion,tandem duplication,inverted tandem duplication,dispersed duplication,dispersed inverted duplication'   # simulated simple SV types
basic_sv_ratio_from_visor = "16:14:14:14:14:14:14"  # ratio of those simple SV types

for chr in chroms:

    working_dir = os.path.join(dir, chr)

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    sv_num_from_visor = int(chr_sv_muns[chr] / (subclone_num))
    # # step 1: get chr's ref
    print("--------------------Start step 1: get chr's ref and index--------------------")

    chr_ref = chr + '.fa'
    chr_ref_path = os.path.join(working_dir, chr_ref)

    cmd_str = "samtools faidx {0} {1} > {2}".format(ref_path, chr, chr_ref_path)
    os.system(cmd_str)
    cmd_str = "samtools faidx {0}".format(chr_ref_path)
    os.system(cmd_str)

    chr_ref_fai_path = os.path.join(working_dir, chr_ref + '.fai')
    chr_snp_h1_path = os.path.join(chr_snp_path, '{0}.h1.fa'.format(chr))

    # step 2: get chr's dim
    print("--------------------Start step 2: get chr's dim--------------------")
    chr_dim_path = os.path.join(working_dir, chr + '.dim.tsv')
    cmd_str = "cut -f1,2 {0} > {1}".format(chr_ref_fai_path, chr_dim_path)
    os.system(cmd_str)

    # step 3: zse VISOR to generate random SV region
    print("--------------------Start step 3: get sv from sim_clone--------------------")
    chr_N_file = os.path.join(N_region_path, 'all_n_region_{0}.txt'.format(chr))
    all_outbed_from_visor = []
    for i in range(subclone_num):
        basic_sv_outbed_from_visor = os.path.join(working_dir, "{0}.sim_clone.all-clone{1}.bed".format(chr, str(i)))
        all_outbed_from_visor.append(basic_sv_outbed_from_visor)

        cmd_str = "Rscript {0}/randomregion.r  -d {1} -n {2} -l {3} -s {4} -v '{5}' -r '{6}' -x {7} > {8}" \
            .format(src_path, chr_dim_path, sv_num_from_visor, sv_len_mean, sv_len_s, basic_sv_type_from_visor,
                    basic_sv_ratio_from_visor, chr_N_file, basic_sv_outbed_from_visor)

        os.system(cmd_str)

    # step 4: use simulate.py add to generate nested sv
    print("--------------------Start step 4: use simulate.py add to generate nested sv--------------------")
    all_nested_bed = []
    for bed in all_outbed_from_visor:
        cmd_str = 'python {0}/simulate.py add -i {1} -w {2} -c {3} -p {4}'.format(src_path, bed, working_dir, chr, chr_snp_h1_path)
        os.system(cmd_str)

        nested_bed = bed.replace('.bed', '') + ".added_SVs.bed"
        all_nested_bed.append(nested_bed)
    all_outbed_from_visor = all_nested_bed
    exit()

    # step 5: use simulate.py sim_clone to add sv to template fa(already has snp)
    all_clone_fa = []
    for i in range(subclone_num):
        sv_bed = all_outbed_from_visor[i]
        clone_fa = chr + '.clone' + str(i)
        all_clone_fa.append(os.path.join(working_dir, clone_fa + '.h1.fa'))

        cmd_str = "python {0}/simulate.py sim_clone -g {1} -bed {2} -o {3} -c {4}" \
            .format(src_path, chr_snp_h1_path, sv_bed, working_dir, clone_fa)
        os.system(cmd_str)

    for fa in all_clone_fa:
        cmd_str = "samtools faidx {0}".format(fa)
        os.system(cmd_str)

    all_clone_dir = []
    for i in range(subclone_num):
        clone_dir = os.path.join(working_dir, 'clone' + str(i))
        all_clone_dir.append(clone_dir)

        cmd_str = "mkdir {0} && mv {1}* {2}".format(clone_dir, all_clone_fa[i], clone_dir)
        os.system(cmd_str)

    print("--------------------Finish step 5: use simulate.py sim_clone to add sv to template fa--------------------")
    
    # # step 6: set purity and bed list
    all_clone_dir_str = ""
    for i in all_clone_dir:
        all_clone_dir_str += i + "/*.fai "
    all_clone_dir_str += chr_ref_fai_path

    purity_out_file = os.path.join(working_dir, 'shorts.laser.multi.bed')

    cmd_str = "cut -f1,2 " + all_clone_dir_str + " | sort | awk '$2 > maxvals[$1] {lines[$1]=$0; maxvals[$1]=$2} END { for (tag in lines) print lines[tag] }' | awk 'OFS=FS=\"\t\"''{print $1, '0', $2, '100.0', '100.0'}' > " + purity_out_file
    os.system(cmd_str)


