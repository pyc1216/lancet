# Fix a problem which causes false negative 


I've been using Lancet to call variants from target  capture data and noticed one **SOMATIC** short tandem repeat be filtered by _HighVafNormal;HighAltCntNormal_  with **WRONG** normal _AD_ which can be detected by Varscan2/Mutect2.

## Before (commit=ce36626c2bb)
1. lancet vcf(fn)

 | CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | normal | tumor | 
 | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
 | chr17 | 37880984 | . | A | ATACGTGATGGCT | 1119.11 | HighVafNormal;HighAltCntNormal | SHARED;FETS=1119.11;TYPE=ins;LEN=12;KMERSIZE=17;SB=14.2495 | GT:AD:SR:SA:DP | 0/1:2101,811:1292,809:2,809:2912 | 0/1:4001,346:2159,1842:181,165:4347 | 

2. varscan2 vcf(tp)

| CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | normal | tumor | 
 | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
| chr17 | 37880984 | . | A | ATACGTGATGGCT | . | PASS | DP=4227;SOMATIC;SS=2;SSC=250;GPV=1E0;SPV=8.049E-26 | GT:GQ:DP:RD:AD:FREQ:DP4 | 0/0:.:1387:1387:0:0%:1068,319,0,0 | 0/1:.:2840:2691:142:5.01%:2287,404,116,26 |


3. mutect2 vcf(tp)

| CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | normal | tumor |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
| chr17 | 37880984 | . | A | ATACGTGATGGCT | . | . | AS_SB_TABLE=2424,2071\|149,149;DP=5352;ECNT=1;MBQ=20,20;MFRL=172,154;MMQ=60,60;MPOS=38;NALOD=3.18;NLOD=348.86;POPAF=6.00;RPA=1,2;RU=TACGTGATGGCT;STR;TLOD=980.44 | GT:AD:AF:DP:F1R2:F2R1:SB | 0/0:1402,0:8.283e-04:1402:713,0:677,0:815,587,0,0 | 0/1:3093,298:0.088:3391:1660,131:1422,134:1609,1484,149,149|

## After

| CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO | FORMAT | normal | tumor |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | 
| chr17 | 37880984 | . | A | ATACGTGATGGCT | 437.116 | PASS | SOMATIC;FETS=437.116;TYPE=ins;LEN=12;KMERSIZE=25;SB=15.3947 | GT:AD:SR:SA:DP | 0/0:1463,0:875,588:0,0:1463 | 0/1:3897,327:2110,1787:168,159:4224 |
