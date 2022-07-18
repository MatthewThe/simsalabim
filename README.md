# SIMSALABIM: Simple Interface for MS Applications

[![PyPI version](https://img.shields.io/pypi/v/simsalabim.svg?logo=pypi&logoColor=FFE873)](https://pypi.org/project/simsalabim/)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/simsalabim.svg?logo=python&logoColor=FFE873)](https://pypi.org/project/simsalabim/)
[![PyPI downloads](https://img.shields.io/pypi/dm/simsalabim.svg)](https://pypistats.org/packages/simsalabim)


## Installation

### With ``pip``

```
pip install simsalabim
```

### From source

```
git clone https://github.com/MatthewThe/simsalabim.git
cd simsalabim
pip install .
```


## Usage

### apl to mzML conversion

The `simsalabim.convert` module automatically detects the output format based on the file extension, e.g. to convert to mzML name the output file `spectra.mzML`.

```
python -m simsalabim.split_apl <mq_andromeda_folder> --output_dir <output_dir>
python -m simsalabim.convert <apl_file> --output_fn <mzml_file>
```

