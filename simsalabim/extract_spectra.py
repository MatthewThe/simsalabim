import os
import sys
from typing import Dict, Optional, List

from pyteomics import mzml


def get_scans(
        file_path: str,
        data: Dict,
        scanidx: Optional[List] = None,
        *args,
        **kwargs
    ) -> None:
    """
    Reads mzml and generates a dataframe containing intensities and m/z values.

    :param file_path: path to a single mzml file.
    :param data: dictionary to be added to by this function
    :param scanidx: optional list of scan numbers to extract. if not specified, all scans will be extracted
    """
    data_iter = mzml.read(source=file_path, *args, **kwargs)
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    for s in scanidx:
        spec = data_iter[f'controllerType=0 controllerNumber=1 scan={s}']
        id = spec['id'].split('scan=')[-1]
        key = f"{file_name}_{id}"
        data[key] = [file_name, id, spec['intensity array'], spec['m/z array']]
    data_iter.close()


if __name__ == "__main__":
    if len(sys.argv) == 3:
        data = {}
        scanidx = list(map(int, sys.argv[2].split(",")))
        get_scans(sys.argv[1], data, scanidx)
        for _, spec_id, intensities, mzs in data.values():
            print(spec_id)
            for m, i in zip(mzs, intensities):
                print(m, i, sep='\t')
    else:
        print("Please specify a mzml file and a comma separated list of scan numbers to extract")
