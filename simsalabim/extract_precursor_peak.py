import sys
import logging

import numpy as np
import matplotlib.pyplot as plt
from biosaur2.utils import MS1OnlyMzML
import mplcursors


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main(argv):
    mzml_file = argv[0]
    
    min_intensity = 0.0
    
    min_mz = float(argv[1])
    max_mz = min_mz + 3.0
    
    min_rt = float(argv[2])
    max_rt = min_rt + 1.0
    
    args = { 'file' : mzml_file,
             'mini' : min_intensity,
             'minmz' : min_mz,
             'maxmz' : max_mz,
             'minrt' : min_rt,
             'maxrt' : max_rt }
    
    logger.info("Extracting")
    
    spectra = process_mzml(args)
    x, y, z = list(), list(), list()
    for spectrum in spectra:
        rt = spectrum['scanList']['scan'][0]['scan start time']
        for mz, i in zip(spectrum['m/z array'], spectrum['intensity array']):
            x.append(rt)
            y.append(mz)
            z.append(np.log10(i))
    
    logger.info("Plotting")
    
    plt.scatter(x, y, c=z)
    plt.xlabel('Retention time [min]')
    plt.ylabel('m/z')
    cbar = plt.colorbar()
    cbar.set_label('log10(intensity)', rotation=90)
    plt.tight_layout()
    
    mplcursors.cursor(hover=True)
    
    plt.show()


def process_mzml(args):

    input_mzml_path = args['file']
    min_intensity = args['mini']
    min_mz = args['minmz']
    max_mz = args['maxmz']
    min_rt = args['minrt']
    max_rt = args['maxrt']

    skipped = 0
    data_for_analyse = []

    cnt = 0

    for z in MS1OnlyMzML(source=input_mzml_path):
        if z['ms level'] == 1:
            if z['scanList']['scan'][0]['scan start time'] < min_rt:
                skipped += 1
                continue
            
            if z['scanList']['scan'][0]['scan start time'] > max_rt:
                skipped += 1
                break

            if 'mean inverse reduced ion mobility array' not in z:
                z['ignore_ion_mobility'] = True
                z['mean inverse reduced ion mobility array'] = np.zeros(len(z['m/z array']))

            idx = z['intensity array'] >= min_intensity
            z['intensity array'] = z['intensity array'][idx]
            z['m/z array'] = z['m/z array'][idx]
            z['mean inverse reduced ion mobility array'] = z['mean inverse reduced ion mobility array'][idx]

            idx = z['m/z array'] >= min_mz
            z['m/z array'] = z['m/z array'][idx]
            z['intensity array'] = z['intensity array'][idx]
            z['mean inverse reduced ion mobility array'] = z['mean inverse reduced ion mobility array'][idx]

            idx = z['m/z array'] <= max_mz
            z['m/z array'] = z['m/z array'][idx]
            z['intensity array'] = z['intensity array'][idx]
            z['mean inverse reduced ion mobility array'] = z['mean inverse reduced ion mobility array'][idx]

            idx = np.argsort(z['m/z array'])
            z['m/z array'] = z['m/z array'][idx]
            z['intensity array'] = z['intensity array'][idx]
            z['mean inverse reduced ion mobility array'] = z['mean inverse reduced ion mobility array'][idx]

            cnt += 1

            # if len(data_for_analyse) > 50:
            #     break

            if len(z['m/z array']):
                data_for_analyse.append(z)
            else:
                skipped += 1


    logger.info('Number of MS1 scans: %d', len(data_for_analyse))
    logger.info('Number of skipped MS1 scans: %d', skipped)

    if len(data_for_analyse) == 0:
        raise Exception('no MS1 scans in input file')

    return data_for_analyse
    

if __name__ == "__main__":
    main(sys.argv[1:])
