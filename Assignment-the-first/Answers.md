# Assignment the First

## Part 1
1. Be sure to upload your Python script. (script name: qscoredist.py)

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
![R1 Qscore Histogram](https://github.com/2020-bgmp/demultiplexing-demiglidden/blob/master/Assignment-the-first/qscore_histRead1_R1.png?raw=true)
![R2 Qscore Histogram](https://github.com/2020-bgmp/demultiplexing-demiglidden/blob/master/Assignment-the-first/qscore_histIndex1_R2.png?raw=true)
![R3 Qscore Histogram](https://github.com/2020-bgmp/demultiplexing-demiglidden/blob/master/Assignment-the-first/qscore_histIndex2_R3.png?raw=true)
![R4 Qscore Histogram](https://github.com/2020-bgmp/demultiplexing-demiglidden/blob/master/Assignment-the-first/qscore_histRead2_R4.png?raw=true)

    2. I suggest a qualtiy score cutoff of 30, which means we would eliminate all data that has a probability of being wrong 1 in a thousand times, or higher. I think that choosing a larger qscore cutoff would be too stringent, but I don't think we should go any lower because we're dealing with over a billion base calls, which can compound error rates fast, so if we kept data with lower qscores, we might end up with a large number of potentially incorrect base calls.
    3. 3,328,051
    command used: '''zcat <file> | sed -n 2~4 | grep "N" | wc -l'''
    
## Part 2
1. Define the problem  
Please see pseudocode_demultiplexing.txt.
2. Describe output  
Please see pseudocode_demultiplexing.txt.
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).  
Please see directories.
4. Pseudocode  
Please see pseudocode_demultiplexing.txt.
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement  
Please see pseudocode_demultiplexing.txt.
