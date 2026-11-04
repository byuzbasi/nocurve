# Raman data and fixed split

The two files come from Cesar Alvarez-Llamas, Vincent Motto-Ros, Christophe
Dujardin, Jérémie Margueritat, David Rodney and David Martin-Calle (2026),
*Raw Raman Spectra of a Heterogeneous Mineral Sample for High-Throughput
Noise2Noise Denoising*, Zenodo v1.
[Source record and DOI 10.5281/zenodo.18244161](https://zenodo.org/records/18244161).
[Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/).

The split uses zero-based source row 440 from `train-data-5ms.h5`. Acquisition
1481 supplies the single training spectrum in `train.tsv`. The channelwise
arithmetic mean of acquisitions 1485–1500 supplies the single test reference
in `test.tsv`. Both files contain all 800 channels in their original order.
The response and reference use source intensity units. The equally spaced
coordinate `x` is used for fitting. The stored `raman_shift` is for plotting.

`split.json` records the source coordinates and file checksums.
