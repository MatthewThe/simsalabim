class Precursor():
    def __init__(self, rtime, precmz, charge, intensity):
        self.precursor = dict()
        self.precursor['mz'] = precmz
        self.precursor['intensity'] = intensity
        self.precursor['charge'] = charge
        self.precursor['rtime'] = rtime # not used

    def precMass(self):
        return helpers.precMassFromPrecMz(self.precursor['mz'], self.precursor['charge'])
