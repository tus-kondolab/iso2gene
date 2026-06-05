# gzip fixtures

`tx2gene_external.tsv.gz` was generated from the adjacent self-authored TSV
with an external gzip implementation, not with iso2gene or miniz:

```bash
gzip -n -9 -c tests/data/gzip/tx2gene_external.tsv > tests/data/gzip/tx2gene_external.tsv.gz
```

SHA256:

```text
21904da6619413572d3cce796dd736583073887d17b6382a1ce0ca5914f5e6d8
```
