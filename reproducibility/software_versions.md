# Software versions

The retained analyses were run with R 4.4.1 and Python 3.11.8 on arm64 macOS. Versions below were taken from the session records generated with the final analysis runs; they are not inferred from current package releases.

## Python

- Python 3.11.8
- h5py (version not captured in the retained session summaries)
- NumPy 2.4.6
- pandas 2.2.3 or 3.0.5, depending on module
- SciPy 1.16.2 or 1.17.1, depending on module

## R and major packages

- R 4.4.1
- data.table 1.18.2
- Matrix 1.7-series
- GSVA 2.0.7
- singscore 1.22.0
- MCPcounter 1.2.0
- xCell 1.1.0
- metafor 4.8-series
- arrow 23.0.1
- hdf5r 1.3.12
- RANN 2.6.2
- survival 3.8-series
- jsonlite 2.0.0
- digest 0.6.39

BayesPrism was used for decomposition in the original analysis environment. Its large fitted objects and private local library are not redistributed. Individual scripts declare their attached packages at the top of each file.
