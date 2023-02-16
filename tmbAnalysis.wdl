version 1.0

workflow tmbAnalysis {
  input {
    File inputMaf
    String intervalFile
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
    description: "Tumour Mutation Burden (TMB) Workflow\n\nTumour Mutation Burden scores are calculated from somatic variant calls coming out of the mutect2 workflow, after annotation with variant effect predictor (vep).  This requires a tumour sample with a matched normal.\n\nThe callable space is the region of the genome where calls are being made.  Accuracy of the TMB metric is dependent on having sufficient depth to ensure coverage over this region. We recommend at a mininum 80X on the tumour and 30X on the matched normal.\n\nTMB is the proportion of the callable space where mutations are identfified.  This is restricted to protein altering mutations of the following vep classes:\n\nMissense_Mutation,In_Frame_Ins,In_Frame_Del,Frame_Shift_Ins,Frame_Shift_Del,Splice_Site,Translation_Start_Site,Nonsense_Mutation,Nonstop_Mutation,Silent"
    dependencies:[]
    output_meta: {
      outputTMB: "output TMB json file"
    }  
  } 
}

task calculateTMB {
  input {
    File inputMaf
    String targetBed
    String outputFileNamePrefix
    String modules = "tmb-r/1.1 python/3.7"
    Int jobMemory = 4
    Int timeout = 6
  }

  parameter_meta {
    inputMaf: "input maf file"
    targetBed: "target bed file"
    outputFileNamePrefix: "prefix for output file"
    modules: "module for running preprocessing"
    jobMemory: "memory allocated to preprocessing, in GB"
    timeout: "timeout in hours"
  }

  command <<<

    set -euo pipefail


    col=$(zcat ~{inputMaf} | grep "Hugo_Symbol" | head -n1 | awk '{ for (i=1; i<=NF; ++i) { if ($i =="FILTER") print i } }')
    
    if [[ $col > 0 ]]; then
      zcat ~{inputMaf} | grep "Hugo_Symbol" | head -n1 > ~{outputFileNamePrefix}.pass.maf
      zcat ~{inputMaf} | awk -F '\t' -v id="$col" '$id=="PASS"' >> ~{outputFileNamePrefix}.pass.maf
    else
      zcat ~{inputMaf} > ~{outputFileNamePrefix}.pass.maf
    fi


    space=$(awk -F '\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM/1000000}' ~{targetBed})
    $TMB_R_ROOT/bin/TMB.R -i ~{outputFileNamePrefix}.pass.maf -o ~{outputFileNamePrefix}.txt -p -c ${space}

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
