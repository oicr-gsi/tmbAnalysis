# tmbAnalysis

workflow to calculate TMB

## Overview

## Dependencies



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
`calculateTMB.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`calculateTMB.timeout`|Int|20|timeout in hours
`calculateTMB.threads`|Int|8|number of cpu threads to be used


### Outputs

Output | Type | Description
---|---|---
`outputTMB`|File|output TMB file


## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
