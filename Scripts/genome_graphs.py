import os
import os.path as p


# Build DBG
# ---------
def build_dbg(interleaved):
    # Move these files into a temp file and pass it to Bifrost
    # For now use CWD for debugging
    cwd = os.getcwd()
    filename = "seqeunce_paths.txt"
    myfile = open(p.join(cwd, filename), 'w')

    myfile.write("\n".join(interleaved))
    myfile.close()


    print("Running Bifrost")

    input_filename = filename
    output_filename = "RSV-bifrost"
    bifrost_str = "Bifrost build \
    -c \
    -k 31 \
    -s {0} \
    -o {1}".format(input_filename, output_filename)
    os.system(bifrost_str)
    return output_filename


# Induce vg
# --------
def induce_vg():
    seqwish_str = "seqwish \
  -s SARS-CoV-2.fasta \
  -p SARS-CoV-2.paf \
  -g SARS-CoV-2.seqwish.gfa"

# Prepare for vg
# --------------

## odgi
def run_odgi():
    odgi_build_str = "odgi build \
  -s \
  -g SARS-CoV-2.seqwish.gfa \
  -o SARS-CoV-2-odgi-graph.vg"

    odgi_chop_str = "odgi chop \
  -i SARS-CoV-2-odgi-graph.vg \
  -c 1024 \
  -o SARS-CoV-2-odgi-chopped.vg"

    odgi_sort_str = "odgi sort \
 -i SARS-CoV-2-odgi-chopped.vg \
 -o SARS-CoV-2-odgi-sorted.vg"

    odgi_view_str = "odgi view \
 -i SARS-CoV-2-odgi-sorted.vg \
 -g \
 > SARS-CoV-2-odgi.gfa"

# vg
# --
def run_vg():
    vg_view_str = "vg view -Fv SARS-CoV-2-odgi.gfa > SARS-CoV-2-vg.vg"

    vg_index_str = "vg index -x SARS-CoV-2.xg -g SARS-CoV-2.gcsa SARS-CoV-2-vg.vg"

    # Map multiple samples against each other
    vg_map_str = "vg map \
    -f ~/projects/Masters/verify/data/reads/simulated/COVID_19/${i}/covid_19_sim_${i}_interleaved.fastq \
    -x SARS-CoV-2.xg \
    -g SARS-CoV-2.gcsa \
    > SARS-CoV-2-${i}.gam"

    vg_coverage = "vg pack \
   -x SARS-CoV-2.xg \
   -g SARS-CoV-2-${i}.gam \
   -d \
   > SARS-CoV-2-${i}.pack.table"


# Analysis
# --------


# Bluntify
 # --------

def bluntify():
    gimbricate_str = "gimbricate -d  \
  -g SARS-CoV-2.gfa \
  -p SARS-CoV-2.paf \
  -f SARS-CoV-2.fasta \
  > SARS-CoV-2.gimbry.gfa"


# Interleave
# ----------
def call_interleave(fp):
    contents = os.listdir(fp)

    if len(contents) == 2:
        f = p.join(fp,contents[0])
        s = p.join(fp,contents[1])
        t = p.join(fp,"interleaved.fq")

        interleave = "interleave-fastq {0} {1} > {2}".format(f, s, t)
        os.system(interleave)

        print("Interleaved %s" % t)
        return t
    else:
        raise Exception("More than 2 files in %s " % fp)

def interleave_reads_in_dir(fp):
    """
    Pass directory with reads.
    """
    interleaved = []
    for d in os.listdir(fp):
       i = call_interleave(p.join(fp,d))
       interleaved.append(i)

    return interleaved

def interleave(fp = os.getcwd()):
    print("Interleaving input files")
    if fp:
        pass
    else:
        fp = os.getcwd()

        # Interleave reads
        return interleave_reads_in_dir(p.join(fp, "data"))

def fetch_reads():
    pass

def main():
    # Read args

    # Interleave reads
    interleaved = interleave()

    # Fetch reads
    fetchReads()

    # Build DBG
    build_dbg(interleaved)

    # Bluntify

    # Induce a variation graph using seqwish

    # Prepare the graph for use with vg

    # vg

main()
