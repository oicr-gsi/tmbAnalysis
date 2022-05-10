# tmbAnalysis

workflow to calculate TMB

## Usage

### Cromwell
```
java -jar cromwell.jar run tmbAnalysis.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputMaf`|File|input maf file
`intervalFile`|String|target bed file
`outputFileNamePrefix`|String|Prefix to use for output file


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`calculateTMB.modules`|String|"tmb-r/1.0"|module for running preprocessing
`calculateTMB.jobMemory`|Int|4|memory allocated to preprocessing, in GB
`calculateTMB.timeout`|Int|6|timeout in hours


### Outputs

Output | Type | Description
---|---|---
`outputTMB`|File|output TMB file


## Commands
 This section lists command(s) run by tmbAnalysis workflow
 
 * Running tmbAnalysis
 
 tmbAnalysis workflow runs the following command (excerpt from .wdl file). inputMaf is a placeholder for an input file.
 
 <<<
 
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
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
