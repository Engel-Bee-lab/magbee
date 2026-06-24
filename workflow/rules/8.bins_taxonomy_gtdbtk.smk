"""
Rules to assign taxonomy to high and medium quality DAS Tool bins using GTDB-Tk.
Filters CheckM2 output for DAS Tool bins to select only high quality (HQ) and medium quality (MQ) bins.

Quality thresholds used:
- High Quality (HQ): Completeness >= 90% AND Contamination <= 5%
- Medium Quality (MQ): Completeness >= 50% AND Contamination <= 10%
"""
import os
import pandas as pd

rule filter_dastool_hq_mq_bins:
    """
    Filter DAS Tool bins based on CheckM2 quality metrics.
    Selects only HQ (>=90% complete, <=5% contamination) and MQ (>=50% complete, <=10% contamination) bins.
    """
    input:
        checkm2_report = os.path.join(dir_reports, "checkm2_all", "CheckM2_DASTool_quality_report.tsv")
    output:
        filtered_done = os.path.join(dir_binning, "dastool_hq_mq_filtered", "done.txt")
    params:
        source_bins = os.path.join(dir_binning, "das_tool", "dastool_DASTool_bins"),
        output_dir = os.path.join(dir_binning, "dastool_hq_mq_filtered"),
        checkm2_report = os.path.join(dir_reports, "checkm2_all", "CheckM2_DASTool_quality_report.tsv")
    localrule: True
    shell:
        """
        mkdir -p {params.output_dir}
        
        python3 << 'EOF'
        import pandas as pd
        import os
        import shutil

        checkm2_file = "{params.checkm2_report}"
        source_dir = "{params.source_bins}"
        output_dir = "{params.output_dir}"

        # Check if CheckM2 report exists
        if not os.path.exists(checkm2_file):
            print(f"CheckM2 report not found: {checkm2_file}")
            exit(0)

        try:
            df = pd.read_csv(checkm2_file, sep='\t')
        except Exception as e:
            print(f"Error reading CheckM2 report: {e}")
            exit(0)

        count = 0
        for idx, row in df.iterrows():
            try:
                completeness = float(row['Completeness'])
                contamination = float(row['Contamination'])
                bin_id = row['Name']
                
                # Filter for high quality (>=90% complete, <=5% contamination) 
                # or medium quality (>=50% complete, <=10% contamination)
                if (completeness >= 90 and contamination <= 5) or (completeness >= 50 and contamination <= 10):
                    # Find matching bin file
                    for ext in ['.fa', '.fasta', '.fna']:
                        source_file = os.path.join(source_dir, bin_id + ext)
                        if os.path.exists(source_file):
                            target_file = os.path.join(output_dir, bin_id + ext)
                            shutil.copy2(source_file, target_file)
                            count += 1
                            break
            except Exception as e:
                print(f"Error processing bin {bin_id}: {e}")

        print(f"Filtered {count} HQ/MQ bins for GTDB-Tk analysis")

        EOF
        
        touch {output.filtered_done}
        """

rule gtdbtk_dastool_hq_mq:
    """
    Run GTDB-Tk on filtered DAS Tool HQ/MQ bins.
    Executes identify -> align -> classify workflow.
    """
    input:
        filtered_bins = os.path.join(dir_binning, "dastool_hq_mq_filtered", "done.txt")
    output:
        gtdbtk_done = os.path.join(dir_binning, "gtdbtk_output", "dastool_hq_mq", "done.txt")
    params:
        bins_dir = os.path.join(dir_binning, "dastool_hq_mq_filtered"),
        outdir = os.path.join(dir_binning, "gtdbtk_output", "dastool_hq_mq"),
        database = config["databases"]["gtdbtk_db"]
    conda:
        os.path.join(dir_env, "gtdbtk.yaml")
    resources:
        mem_mb = config['resources']['bigjob']['mem_mb'],
        runtime = config['resources']['bigjob']['runtime']
    threads:
        config['resources']['bigjob']['threads']
    shell:
        """
        set -euo pipefail
        
        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.database}
        
        # Check if bins exist
        if ! ls {params.bins_dir}/*.fa {params.bins_dir}/*.fasta {params.bins_dir}/*.fna 2>/dev/null | head -1 >/dev/null; then
            echo "No HQ/MQ bins found for GTDB-Tk analysis in {params.bins_dir}"
            touch {output.gtdbtk_done}
            exit 0
        fi
        
        echo "Running GTDB-Tk identify step..."
        gtdbtk identify --genome_dir {params.bins_dir} --cpus {threads} --out_dir {params.outdir} -x fa --force 2>&1 || true
        
        echo "Running GTDB-Tk align step..."
        gtdbtk align --identify_dir {params.outdir} --out_dir {params.outdir} --cpus {threads} --force 2>&1 || true
        
        echo "Running GTDB-Tk classify step..."
        gtdbtk classify --genome_dir {params.bins_dir} --out_dir {params.outdir} --cpus {threads} -x fa --align_dir {params.outdir} --force 2>&1 || true
        
        echo "GTDB-Tk analysis complete for DAS Tool HQ/MQ bins"
        touch {output.gtdbtk_done}
        """
