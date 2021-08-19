## Minimal cDNA-detector batching wrapper

**cDNA-detector**: [preprint](https://www.biorxiv.org/content/10.1101/2021.08.11.455962v1.full) | [GitHub](https://github.com/rheinbaylab/cDNA-detector)

### Setup:

Install [cDNA-detector](https://github.com/rheinbaylab/cDNA-detector) from the [Rheinbay Lab](https://www.massgeneral.org/cancer-center/clinical-trials-and-research/center-for-cancer-research/investigators/rheinbay-lab).

I've cloned the cDNA-detector repository adjacent to repository from which I am executing the commands in this wrapper. Thus, this wrapper references the main python script in that library:
```
/home/mvinyard/software/cDNA-detector/cdna-detector.py
```


### Usage:


**Define the path to all 10x samples**:

Assumes the following structure where only the path before what is shown in brackets is passed to the module:
```python=
data_path = "/home/mvinyard/data/10x_samples/[SAMPLE/outs/possorted_bam.bam]"
```

**Define output path**:

Only a single directory (or none if the working directory suits you) need be definied. Sample names do not need to be defined; this is done by the wrapper class when each sample is detected and outputs are stored separately. 
```python=
outdir = "/path/to/cDNA_detector/outs/[SAMPLE]/"
```


```python=
cDNA = cDNA_detector()
cDNA.preflight(data_path=data_path, outdir=outdir, gene_model="hg38")
```

After running the preflight setup, one can check the in and out paths as well as IDs for each sample:
```python=
cDNA.SampleDict()
```
```
{'Sample_1': ['/path/to/data/Sample_1/outs/possorted_bam.bam',
  '/home/mvinyard/results/cDNA-detector/Sample_1'],
 'Sample_2': ['/path/to/data/Sample_2/outs/possorted_bam.bam',
  '/home/mvinyard/results/cDNA-detector/Sample_2'],
 'Sample_3': ['/path/to/data/Sample_3/outs/possorted_bam.bam',
  '/home/mvinyard/results/cDNA-detector/Sample_3'],
 'Sample_4': ['/path/to/data/Sample_4/outs/possorted_bam.bam',
  '/home/mvinyard/results/cDNA-detector/Sample_4'],
              .
              .
              .
 'Sample_N': ['/path/to/data/Sample_N/outs/possorted_bam.bam',
  '/home/mvinyard/results/cDNA-detector/Sample_N']}

```

**Finally, Run the `cDNA-detector` Module**:
```
cDNA.detect()
```


### Notes:
- other modules of cDNA-detector, namely `prepare` and `clean` have not yet been implemented in this wrapper. 
