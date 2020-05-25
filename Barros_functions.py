import numpy as np
import physconst as pc
import norm

# NOTE: Gravitational constant is in units of kpc, Msun, Myr
# The values in the table files are in units of 10^10 Msun and kpc

nm = norm.PhysNorm(m=pc.msun, x=pc.kpc, t=pc.Myr, temp=1., curr=1.)


def miyamoto_nagai_density_m1(r, z, m, a, b):

    zeta = np.sqrt(z * z + b * b)
    apz2 = (a + zeta) * (a + zeta)

    return ( (b * b * m / (4 * np.pi)) * (a * r * r + (a + 3 * zeta) * (a + zeta) * (a + zeta)) /
             (zeta * zeta * zeta * (r * r + apz2) ** 2.5) )


def miyamoto_nagai_density_m2(r, z, m, a, b):

    zeta = np.sqrt(z * z + b * b)
    apz2 = (a + zeta) * (a + zeta)

    return ( (b * b * m / (4 * np.pi)) * 3 * (a + zeta) *
             (r * r * (zeta * zeta - a * zeta + a * a) + apz2 * (zeta * zeta + 4 * a * zeta + a * a)) /
             (zeta * zeta * zeta * (r * r + apz2) ** 3.5) )


def miyamoto_nagai_density_m3(r, z, m, a, b):

    zeta = np.sqrt(z * z + b * b)
    apz2 = (a + zeta) * (a + zeta)
    zeta3 = zeta * zeta * zeta

    return ( (b * b * m / (4 * np.pi)) * (
            3 * r * r * r * r * zeta3 +
            r * r * apz2 * (6 * zeta3  + 15 * a * zeta * zeta - 10 * a * a * zeta + 5 * a * a * a) +
            apz2 * apz2 * (3 * zeta3 + 15 * a * zeta * zeta + 25 * a * a * zeta + 5 * a * a * a) ) /
            (zeta * zeta * zeta * (r * r + apz2) ** 4.5) )


def triple_miyamoto_nagai_density_m1(r, z, m1, a1, m2, a2, m3, a3, b):

    return miyamoto_nagai_density_m1(r, z, m1, a1, b) + \
           miyamoto_nagai_density_m1(r, z, m2, a2, b) + \
           miyamoto_nagai_density_m1(r, z, m3, a3, b)


def triple_miyamoto_nagai_density_m2(r, z, m1, a1, m2, a2, m3, a3, b):

    return miyamoto_nagai_density_m2(r, z, m1, a1, b) + \
           miyamoto_nagai_density_m2(r, z, m2, a2, b) + \
           miyamoto_nagai_density_m2(r, z, m3, a3, b)


def triple_miyamoto_nagai_density_m3(r, z, m1, a1, m2, a2, m3, a3, b):

    return miyamoto_nagai_density_m3(r, z, m1, a1, b) + \
           miyamoto_nagai_density_m3(r, z, m2, a2, b) + \
           miyamoto_nagai_density_m3(r, z, m3, a3, b)


def miyamoto_nagai_potential_m1(r, z, m, a, b):

    zeta = np.sqrt(z * z + b * b)
    apz2 = (a + zeta) * (a + zeta)

    Gnm = pc.newton / nm.newton

    return -Gnm * m / np.sqrt(r * r + apz2)


def miyamoto_nagai_potential_m2(r, z, m, a, b):

    zeta = np.sqrt(z * z + b * b)
    apz2 = (a + zeta) * (a + zeta)

    return miyamoto_nagai_potential_m1(r, z, m, a, b) * (1 + a * (a + zeta) / (r * r + apz2))


def miyamoto_nagai_potential_m3(r, z, m, a, b):

    zeta = np.sqrt(z * z + b * b)
    apz2 = (a + zeta) * (a + zeta)

    return ( miyamoto_nagai_potential_m2(r, z, m, a, b) * (1 + a * (a + zeta) / (r * r + apz2)) -
             miyamoto_nagai_potential_m1(r, z, m, a, b) * 0.333333 * a * a * (r * r - 2 * apz2) /
             ((r * r  + apz2) * (r * r  + apz2)) )


def triple_miyamoto_nagai_potential_m1(r, z, m1, a1, m2, a2, m3, a3, b):

    return miyamoto_nagai_potential_m1(r, z, m1, a1, b) + \
           miyamoto_nagai_potential_m1(r, z, m2, a2, b) + \
           miyamoto_nagai_potential_m1(r, z, m3, a3, b)


def triple_miyamoto_nagai_potential_m2(r, z, m1, a1, m2, a2, m3, a3, b):

    return miyamoto_nagai_potential_m2(r, z, m1, a1, b) + \
           miyamoto_nagai_potential_m2(r, z, m2, a2, b) + \
           miyamoto_nagai_potential_m2(r, z, m3, a3, b)


def triple_miyamoto_nagai_potential_m3(r, z, m1, a1, m2, a2, m3, a3, b):

    return miyamoto_nagai_potential_m3(r, z, m1, a1, b) + \
           miyamoto_nagai_potential_m3(r, z, m2, a2, b) + \
           miyamoto_nagai_potential_m3(r, z, m3, a3, b)


def hernquist_density(r, z, m, a):

    rr = np.sqrt(r * r + z * z)

    return m / (2. * np.pi) * a / (rr * (rr + a) * (rr + a) * (rr + a))


def hernquist_potential(r, z, m, a):

    rr = np.sqrt(r * r + z * z)
    Gnm = pc.newton / nm.newton

    return -Gnm * m / (rr + a)


def logarithmic_density(r, z, rh, vh):

    Gnm = pc.newton / nm.newton

    return ( vh * vh / (4. * np.pi * Gnm) * (r * r + z * z + 3 * rh * rh) /
             ((r * r + z * z + rh * rh) * (r * r + z * z + rh * rh)) )


def logarithmic_potential(r, z, rh, vh):

    return 0.5 * vh * vh * np.log(r * r + z * z + rh * rh)


class BarrosMW:

    def __init__(self, file):
        """
        Parameters
        ----------
        file   File to read in

        Read and store the parameters.
        Only for Model MI in their paper. (I.e., not the one with the ring structure)

        Store data in units of kpc and Msun. Masses in table must be converted from 10^10 Msun to 1 Msun.
        """

        # This is the order of the parameters in the pot potential files
        self.disc_params = ['M1', 'a1', 'M2', 'a2', 'M3', 'a3', 'b']
        self.sph_params = ['Mb', 'ab']
        self.DM_params = ['rh', 'vh']

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

                # Convert masses from units of 10^10 Msun to 1 Msun
                for i_param in range(len(self.disc_params) - 2):
                    if i_param % 2 == 0:
                        p_discs[i_disc][i_param] *= 1.e10

            # Collect spheroids
            n_sph = int(fh.readline())
            p_sph = []
            for i_sph in range(n_sph):
                p_sph += [[float(n) for n in fh.readline().split(' ')]]
            p_sph[0][0] *= 1.e10       # into Msun
            p_sph[1][1] *= 1.e5 / nm.v  # into kpc / Myr

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

            # Units converter
            self.kpc_Myr_to_km_s = 977.77

    def thin_disc_density(self, r, z):
        return triple_miyamoto_nagai_density_m3(r, z, *self.p_thin)

    def thin_disc_potential(self, r, z):
        return triple_miyamoto_nagai_potential_m3(r, z, *self.p_thin)

    def thick_disc_density(self, r, z):
        return triple_miyamoto_nagai_density_m1(r, z, *self.p_thick)

    def thick_disc_potential(self, r, z):
        return triple_miyamoto_nagai_potential_m1(r, z, *self.p_thick)

    def HI_disc_density(self, r, z):
        return triple_miyamoto_nagai_density_m2(r, z, *self.p_HI)

    def HI_disc_potential(self, r, z):
        return triple_miyamoto_nagai_potential_m2(r, z, *self.p_HI)

    def H2_disc_density(self, r, z):
        return triple_miyamoto_nagai_density_m3(r, z, *self.p_H2)

    def H2_disc_potential(self, r, z):
        return triple_miyamoto_nagai_potential_m3(r, z, *self.p_H2)

    def bulge_sph_density(self, r, z):
        return hernquist_density(r, z, *self.p_bulge)

    def bulge_sph_potential(self, r, z):
        return hernquist_potential(r, z, *self.p_bulge)

    def DM_sph_density(self, r, z):
        return logarithmic_density(r, z, *self.p_DM)

    def DM_sph_potential(self, r, z):
        return logarithmic_potential(r, z, *self.p_DM)

    def total_potential(self, r, z):
        return self.thin_disc_potential(r, z) + self.thick_disc_potential(r, z) + \
               self.HI_disc_potential(r, z) + self.H2_disc_potential(r, z) + \
               self.bulge_sph_potential(r, z) + self.DM_sph_potential(r, z)

    def total_density(self, r, z):
        return self.thin_disc_density(r, z) + self.thick_disc_density(r, z) + \
               self.HI_disc_density(r, z) + self.H2_disc_density(r, z) + \
               self.bulge_sph_density(r, z) + self.DM_sph_density(r, z)
