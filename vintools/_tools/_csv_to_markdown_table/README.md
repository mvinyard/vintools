# CSV to markdown table

#### .csv file input:
```
sample,index,metadata
sample_hcz_105,SI-NA-D1,control
sample_hcz_106,SI-NA-A2,control
sample_hcz_107,SI-NA-D2,control
sample_hcz_108,SI-NA-B2,control
sample_hcz_109,SI-NA-C3,exp
sample_hcz_110,SI-NA-B3,exp
sample_hcz_111,SI-NA-C3,exp
sample_hcz_112,SI-NA-D3,exp
```

#### python utility

There are two mdoes of making this conversion using `vintools`:

1. Instantiate a class that allows you to view the formatted table. Most useful for when working in a notebook:
    
    
    ```python
    import vintools as v
    
    MarkdownTable = v.tl.create_md_table()
    ```
    
    ```python
    MarkdownTable.read_df(infile)
    MarkdownTable.df.head()
    ```
    
    ```python
    MarkdownTable.make_table()
    MarkdownTable.view()
    ```
    
    >|sample|index|metadata|
    |--|--|--|
    |sample_hcz_105|SI-NA-D1|control|
    |sample_hcz_106|SI-NA-A2|control|
    |sample_hcz_107|SI-NA-D2|control|
    |sample_hcz_108|SI-NA-B2|control|
    |sample_hcz_109|SI-NA-C3|exp|
    |sample_hcz_110|SI-NA-B3|exp|
    |sample_hcz_111|SI-NA-C3|exp|
    |sample_hcz_112|SI-NA-D3|exp|
    
    
    ```python
    MarkdownTable.write(outfile)
    ```
    ><input>.csv converted to markdown table and saved: <output>.md


2. Run a one-line function: 

    ```python
    v.tl.quicktable(infile, outfile="table_quick.md")
    ```
    ><input>.csv converted to markdown table and saved: <output>.md


* [**Example python notebook**](https://github.com/mvinyard/vintools/blob/main/notebooks/csv_to_markdown_table.ipynb)


#### markdown-beautified output:

|sample|index|metadata|
|--|--|--|
|sample_hcz_105|SI-NA-D1|control|
|sample_hcz_106|SI-NA-A2|control|
|sample_hcz_107|SI-NA-D2|control|
|sample_hcz_108|SI-NA-B2|control|
|sample_hcz_109|SI-NA-C3|exp|
|sample_hcz_110|SI-NA-B3|exp|
|sample_hcz_111|SI-NA-C3|exp|
|sample_hcz_112|SI-NA-D3|exp|
