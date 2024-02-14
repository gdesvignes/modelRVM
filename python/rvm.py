import numpy as np

def rvm(self, x, zeta, beta, phi0, psi0):
        """Compute the RVM curve given the geometrical angles

        """
        x = np.deg2rad(x)
        zeta = np.deg2rad(zeta)
        beta = np.deg2rad(beta)
        phi0 = np.deg2rad(phi0)
        psi0 = np.deg2rad(psi0)

        argx = np.cos(zeta-beta)*np.sin(zeta) - np.sin(zeta-beta)*np.cos(zeta)*np.cos(x-phi0)
        argy =  np.sin(zeta-beta) * np.sin(x-phi0)
        return -np.arctan(argy/argx) + psi0
