# vintools utility functions


## Data manipulation functions

### Parse genomics datafiles

* separate interleafed fastq files:

    ```python
    v.ut.deinterleaf(interleafed_file, outfile_prefix="./",)
    ```

* read a **`.vcf`** file: 
    ```python
    v.ut.read_vcf(path)
    ```
* parse a GTF by genes:
    ```python
    v.ut.parse_gtf(gtf_file_path)
    ```
* Interacting with devices and pytorch data structures:
    ```python
    v.ut.set_device(gpu="0")
    v.ut.torch_device(array, device=_set_device())
    ```
    
### Other data manipulation functions: 
* python-gsutils interface
* flatten a list
* get data bounds
* smooth data points


## System-interaction functions

```python
v.ut.imports()
```

```python
v.ut.version()
```

```python
v.ut.fetch_n_cpus()
```

```python
v.ut.use_n_cores()
```

```python
v.ut.git_clone()
```

```python
v.ut.which()
```

```python
v.ut.mkdir_flex()
```

```python
v.ut.findNdestroy_pycache()
```

```python
v.ut.email()
```

```python
v.ut.import_from_string()
```

### UX functions
```python
v.ut.format_pystring()
```

```python
v.ut.done()
```