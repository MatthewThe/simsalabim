# SIMSALABIM: Simple Interface for MS Applications

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

