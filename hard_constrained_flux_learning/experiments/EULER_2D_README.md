2D Euler source is stored as a gzip-compressed Python file because the connector
used in this research session has a single-file text payload limit.

After cloning:

```bash
gunzip -k experiments/euler_2d_hcfl.py.gz
python experiments/euler_2d_hcfl.py --seed 0
```

The decompressed file is the exact experiment source used for the reported
three-seed results.
