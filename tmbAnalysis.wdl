version 1.0

workflow tmbAnalysis {
  input {
    File inputMaf
    String? intervalFile
    String outputFileNamePrefix
  }

  call calculateTMB {
    input: 
      inputMaf = inputMaf,
      targetBed = intervalFile,
      outputFileNamePrefix = outputFileNamePrefix
  }

  output {
    File outputTMB = calculateTMB.outputTMB
  }  

  parameter_meta {
    inputMaf: "input maf file"
    intervalFile: "target bed file"
    outputFileNamePrefix: "Prefix to use for output file"
  }

  meta {
    author: "Xuemei Luo"
    email: "xuemei.luo@oicr.on.ca"
    description: "Tumour Mutation Burden (TMB) Workflow\n\nTumour Mutation Burden scores are calculated from somatic variant calls coming out of the mutect2 workflow, after annotation with variant effect predictor (vep).  This requires a tumour sample with a matched normal.\n\nAccuracy of the TMB metric is dependent on having sufficient depth to ensure coverage over this region. We recommend at a mininum 80X on the tumour and 30X on the matched normal.\n\nThe TMB is the number of PASS somatic mutations with the VAF > 10% per megabase within the target region.\n\nThe exome mutations are the mutations within the exome region, and the mutations are restricted to protein altering mutations of the following vep classes:\n\nMissense_Mutation,In_Frame_Ins,In_Frame_Del,Frame_Shift_Ins,Frame_Shift_Del,Splice_Site,Translation_Start_Site,Nonsense_Mutation,Nonstop_Mutation.\n\nThe Genome mutations are all mutations within the genome region."
    dependencies:[]
    output_meta: {
    outputTMB: {
        description: "output TMB json file",
        vidarr_label: "outputTMB"
    }
}
  } 
}

task calculateTMB {
  input {
    File inputMaf
    String? targetBed
    Float exome_targetSpace = 34.0
    Float genome_targetSpace = 3095.978588
    Float minVAF = 0.1
    Boolean doGenomicTMB = false
    String outputFileNamePrefix
    String modules = "tmb-r/1.2 python/3.7"
    Int jobMemory = 4
    Int timeout = 6
  }

  parameter_meta {
    inputMaf: "input maf file"
    targetBed: "target bed file."
    exome_targetSpace: "exome target region in Megabase. It is used to calculate exome TMB only when targetBed is not provided"
    genome_targetSpace: "genome region in Megabase. It is used to calculate Genomic TMB"
    outputFileNamePrefix: "prefix for output file"
    minVAF: "minimum VAF value"
    doGenomicTMB: "if it is true, calculate Genomic TMB"
    modules: "module for running preprocessing"
    jobMemory: "memory allocated to preprocessing, in GB"
    timeout: "timeout in hours"
  }

  command <<<

    set -euo pipefail


    col_pass=$(zcat ~{inputMaf} | grep "Hugo_Symbol" | head -n1 | awk '{ for (i=1; i<=NF; ++i) { if ($i =="FILTER") print i } }')
    col_t_depth=$(zcat ~{inputMaf} | grep "Hugo_Symbol" | head -n1 | awk '{ for (i=1; i<=NF; ++i) { if ($i =="t_depth") print i } }')
    col_t_alt_count=$(zcat ~{inputMaf} | grep "Hugo_Symbol" | head -n1 | awk '{ for (i=1; i<=NF; ++i) { if ($i =="t_alt_count") print i } }')

    ## filter the mutations with PASS and minVAF
    if [[ $col_pass > 0  && col_t_depth > 0 && col_t_alt_count > 0 ]]; then
      zcat ~{inputMaf} | grep "Hugo_Symbol" | head -n1 > ~{outputFileNamePrefix}.pass.maf
      zcat ~{inputMaf} | awk -F '\t' -v pass_val="$col_pass" '$pass_val=="PASS"' | awk -F '\t' -v t_depth="$col_t_depth" -v alt_count="$col_t_alt_count" -v minVaf=~{minVAF} '$alt_count > 0 && $t_depth > 0 && $alt_count/$t_depth >= minVaf'  >> ~{outputFileNamePrefix}.pass.maf
    else
      zcat ~{inputMaf} > ~{outputFileNamePrefix}.pass.maf
    fi

    ## calculate target space if a target bed file is provided
    if [ -z "~{targetBed}" ]; then
        space=~{exome_targetSpace};
    else space=$(awk -F '\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM/1000000}' ~{targetBed})
    fi

    ## calculate exome TMB, subset to protein altering mutations
    $TMB_R_ROOT/bin/TMB.R -i ~{outputFileNamePrefix}.pass.maf -o ~{outputFileNamePrefix}.tmp1.txt -p -c ${space}

    cat ~{outputFileNamePrefix}.tmp1.txt | sed 's/Total_Mutations/Exome_Total_Mutations/' | sed 's/Mutation_burden/Exome_Mutation_burden/' | sed 's/Callable_space/Exome_Callable_space/' > ~{outputFileNamePrefix}.exome.txt

    ## calculate genome TBM if doGenomicTMB is true
    if ~{doGenomicTMB} ; then
      $TMB_R_ROOT/bin/TMB.R -i ~{outputFileNamePrefix}.pass.maf -o ~{outputFileNamePrefix}.tmp2.txt -c ~{genome_targetSpace}
      cat ~{outputFileNamePrefix}.tmp2.txt | sed 's/Total_Mutations/Genome_Total_Mutations/' | sed 's/Mutation_burden/Genome_Mutation_burden/' | sed 's/Callable_space/Genome_Callable_space/' > ~{outputFileNamePrefix}.genome.txt
      join -t $'\t' ~{outputFileNamePrefix}.exome.txt ~{outputFileNamePrefix}.genome.txt > ~{outputFileNamePrefix}.txt
    else
      cat ~{outputFileNamePrefix}.exome.txt > ~{outputFileNamePrefix}.txt
    fi

    python3 <<CODE
    import json

    metrics = []
    file = open('~{outputFileNamePrefix}.txt', 'r')
    a = file.readline()
    titles = [t.strip() for t in a.split('\t')]
    for line in file:
        d = {}
        for t, f in zip(titles, line.split('\t')):
            d[t] = f.strip()
        metrics.append(d)
    with open('~{outputFileNamePrefix}.json', 'w') as json_file:
        json.dump(metrics, json_file)

    CODE

  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File outputTMB = "~{outputFileNamePrefix}.json"
  }

}
