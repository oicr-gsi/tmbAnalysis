version 1.0

workflow tmbAnalysis {
  input {
    File inputMaf
    String intervalFile
    String outputFileNamePrefix
}

  parameter_meta {
    inputMaf: "input maf file"
    intervalFile: "target bed file"
    outputFileNamePrefix: "Prefix to use for output file"
  }
   
  call calculateTMB {
    input: 
      inputMaf = inputMaf,
      targetBed = intervalFile,
      outputFileNamePrefix = outputFileNamePrefix
  }

  meta {
    author: "Xuemei Luo"
    email: "xuemei.luo@oicr.on.ca"
    description: "workflow to calculate TMB"
    dependencies:
    []
  }
  
  output {
    File outputTMB = calculateTMB.outputTMB
  }
}

task calculateTMB {
  input {
   File inputMaf
   String targetBed
   String outputFileNamePrefix
   String modules = "tmb-r/1.0"
   Int jobMemory = 24
   Int timeout = 20
   Int threads = 8
  }

  parameter_meta {
  inputMaf: "input maf file"
  targetBed: "target bed file"
  outputFileNamePrefix: "prefix for output file"
  modules: "module for running preprocessing"
  jobMemory: "memory allocated to preprocessing, in GB"
  timeout: "timeout in hours"
  threads: "number of cpu threads to be used"
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

  >>>

  runtime {
    memory:  "~{jobMemory} GB"
    modules: "~{modules}"
    cpu:     "~{threads}"
    timeout: "~{timeout}"
  }

  output {
    File outputTMB = "~{outputFileNamePrefix}.txt"
  }

}
