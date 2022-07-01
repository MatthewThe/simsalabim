class Spectrum:
    def __init__(self, scannr, precursors, fragment_ions, raw_file):
        self.spectrum = dict()
        self.spectrum["precursor_information"] = precursors
        self.spectrum["id"] = f"scan={scannr}"
        self.spectrum["params"] = [{'ms level': 2}]
        
        self.spectrum["mz_array"], self.spectrum["intensity_array"] = zip(*fragment_ions)
        
        self.raw_file = raw_file
    
    def get_ms_level(self):
        for d in self.spectrum["params"]:
            for k, v in d.items():
                if k == "ms level":
                    return v
        return -1
