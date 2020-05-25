import numpy as np


def power_law_spheroid(r, z, rho0, q, gamma, beta, r0, rcut):
    rp = np.sqrt(r * r + (z * z / (q * q)))
    return rho0 / ((rp / r0) ** gamma + (1. + (rp / r0) ** (beta - gamma))) * np.exp(-(rp / rcut) ** 2)


def exponential_disc(r, z, sigma0, rd, zd, r0, eps):

    # z component, dependent on value of zd
    if zd < 0:
        dz = 1. / (np.cosh(z / (2. * zd)) ** 2 * (-4. * zd))
    else:
        dz = np.exp(- np.fabs(z) / zd) / (2. * zd)

    return sigma0 * np.exp(-r0 / r - r / rd + eps * np.cos(np.pi * r / rd)) * dz


class PJM_MW:

    def __init__(self, file):
        """
        Parameters
        ----------
        file   File to read in

        Read and store the parameters
        """

        # This is the order of the parameters in the Tpot potential files
        self.disc_params = ['Sigma0', 'Rd', 'zd', 'R0', 'eps']
        self.sph_params = ['rho0', 'q', 'gamma', 'beta', 'r0', 'rcut']

        # This is the order of the components in the Tpot potential files
        self.disc_labels = ['thin', 'thick', 'HI', 'H2']
        self.sph_labels = ['bulge', 'DM']

        # Read file and store values
        with open(file, 'r') as fh:

            # Collect discs
            n_discs = int(fh.readline())
            p_discs = []
            for i_disc in range(n_discs):
                p_discs += [[float(n) for n in fh.readline().split(' ')]]

            # Collect spheroids
            n_sph = int(fh.readline())
            p_sph = []
            for i_sph in range(n_sph):
                p_sph += [[float(n) for n in fh.readline().split(' ')]]

            # Explicitly store values for the components
            self.p_thin = p_discs[0]
            self.p_thick = p_discs[1]
            self.p_HI = p_discs[2]
            self.p_H2 = p_discs[3]
            self.p_bulge = p_sph[0]
            self.p_DM = p_sph[1]

            # Store other parameters as well
            self.n_discs = n_discs
            self.p_discs = p_discs
            self.n_sph = n_sph
            self.p_sph = p_sph

    def thin_disc(self, r, z):
        return exponential_disc(r, z, *self.p_thin)

    def thick_disc(self, r, z):
        return exponential_disc(r, z, *self.p_thick)

    def HI_disc(self, r, z):
        return exponential_disc(r, z, *self.p_HI)

    def H2_disc(self, r, z):
        return exponential_disc(r, z, *self.p_H2)

    def bulge_sph(self, r, z):
        return power_law_spheroid(r, z, *self.p_bulge)

    def DM_sph(self, r, z):
        return power_law_spheroid(r, z, *self.p_DM)

