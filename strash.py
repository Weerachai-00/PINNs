import scipy.integrate as integrate
import numpy as np
import numpy.random as rd
from fractions import Fraction
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import scipy as sp
import time
import configparser
from sys import argv
import diagonalisation
import eoms
import output
import plot
import sys
from sys import argv
import model_class
import matplotlib.pyplot as plt

class hubble_calculator(object):

    mpl = 2.435e27

    def __init__(self, xin):
        """Initialize the object with parameter array xin, which can alternatively be specified as a string which points to the file name of a configuration.ini file.  Samples initial conditions and masses via the diagonalisation routine."""

        if isinstance(xin, str):              # มีการแก้ไขในส่วนตรงนี้ใหม่ ให้มีการตรวจสอบว่า xin เป็นข้อมูลประเภท string จริงหรือไม่
            xin = self.read_configfile(xin)

        self.unpack(xin)

        self.ma_array, self.fef, self.phiin_array, self.phidotin_array = diagonalisation.diag(
            self.mo, self.c, self.n, self.a0, self.b0, self.sa, self.sb, self.s1, self.s2, self.phi_range, self.phidotin)

        self.rhoin_array = eoms.rhoinitial(self.phidotin_array, self.phiin_array, self.ma_array, self.n)
        self.y0 = eoms.yinitial(self.n, self.phiin_array, self.phidotin_array, self.rhoin_array, self.ain)

    def unpack(self, xin):
        """Simple routine to take an array of parameters xin and maps it to a series of class attributes.  Defines appropriate ordering between attributes and a vector of inputs."""

        (self.mo, self.c, self.n, self.N, self.n_cross, self.a0, self.sa, self.b0, self.sb, self.s1, self.s2, 
         self.phi_range, self.phidotin, self.ain, self.tin, self.tfi, self.rhocrit, self.rho_m0, self.rho_r0, self.rhol) = xin

    def eq(self, y, t):
        """Equations of motion."""
        return eoms.deriv_wfromphi(y, t, self.n, self.n_cross, self.ma_array, self.rho_m0, self.rho_r0, self.rhol)

    def solver(self):
        """Solve the equations of motion with initial and final time set by class attributes."""
        self.t = np.linspace(self.tin, self.tfi, self.N)
        self.y = integrate.odeint(self.eq, self.y0, self.t, mxstep=1000000000)
        print(self.y)
        uuut = self.t 
        uuuy = self.y
        return uuut, uuuy
    def output(self):
        """Obtain some derived quantities."""
        N = len(self.t)
        self.rhoa = output.axionrho(self.y, N, self.n)
        self.P = output.pressure(self.y, self.ma_array, N, self.n)
        self.w = output.w(self.P, self.rhoa, N)
        self.rhom, self.rhor = output.dense(self.rho_m0, self.rho_r0, N, self.y)
        self.rhosum = output.totalrho(self.rhom, self.rhol, self.rhor, self.rhoa, N)
        self.H = output.hubble(self.t, self.rhosum)
        self.z = output.redshift(self.y, N)
        self.ODM, self.ODE = output.dark(self.rhom, self.rhor, self.rhol, self.rhoa, self.ma_array, self.rhosum, self.n, self.y)

    def printout(self):
        """Print calculated quantities to stdout."""
        print("H", "w", "z")
        for i in range(len(self.H)):
            print(self.H[i], self.z[i], self.w[i])
        if __name__ == "__main__":
            if len(sys.argv) > 2 and sys.argv[2] == "cosmo":             # แก้ไขตรงนี้เพิ่มเติมเพื่อให้ทำงานง่ายขึ้นโดยเปลี่ยนเป็น 1 if  และแก้ไขตรง len(sys.argv) > 1 เป็น len(sys.argv) > 2 
                                                                         # การตรวจสอบนี้เพื่อดูว่ามีอาร์กิวเมนต์ที่สอง (sys.argv[2]) ถูกส่งเข้ามาหรือไม่ และตรวจสอบว่ามันมีค่าตรงกับ "cosmo" และยังเป็ฯการเช็คอีกว่า argv > 2 จริงๆ
                plot.cosmology(self.rhoa, self.rhosum, self.rhor, self.rhom, self.y)

    def read_configfile(self, fname):
        """Read a configuration file to set parameters for the calculation."""

        config = configparser.RawConfigParser()
        config.read(fname)             # แทนที่จะให้อ่าน argv[1] เปลี่ยนเป็น fname เพื่อความเข้าใจของตัวเอง

        mo = config.getint('Model_Selection', 'Model')
        c = config.getfloat('Model_Selection', 'Dimension')
        n = config.getint('General', 'Number of Axions')
        N = config.getint('General', 'Number of time steps')
        n_cross = config.getint('General', 'Number of Crossings')

        a0 = config.getfloat('Hyperparameter', 'a0')
        sa = config.getfloat('Hyperparameter', 'sigma_a')
        b0 = config.getfloat('Hyperparameter', 'b0')
        sb = config.getfloat('Hyperparameter', 'sigma_b')
        s1 = config.getfloat('Hyperparameter', 's1')
        s2 = config.getfloat('Hyperparameter', 's2')

        phi_range = config.getfloat('Initial Conditions', 'phi_in_range')
        phidotin = config.getfloat('Initial Conditions', 'phidot_in')
        ain = config.getfloat('Initial Conditions', 'a_in')
        tin = config.getfloat('Initial Conditions', 't_in')
        tfi = config.getfloat('Initial Conditions', 't_fi')

        rhocrit = config.getfloat('Physical Constants', 'rho_critical')
        rho_m0 = config.getfloat('Physical Constants', 'rho_matter_0')
        rho_r0 = config.getfloat('Physical Constants', 'rho_radiation_0')
        rhol = config.getfloat('Physical Constants', 'rho_lambda')

        xin = [mo, c, n, N, n_cross, a0, sa, b0, sb, s1, s2, phi_range, phidotin, ain, tin, tfi, rhocrit, rho_m0, rho_r0, rhol]
        print("mo, c, n, N, n_cross, a0, sa, b0, sb, s1, s2, phi_range, phidotin, ain, tin, tfi, rhocrit, rho_m0, rho_r0, rhol")
        print(xin)
        return xin


def main():
    if len(argv) < 2:           
        raise Exception('Need to specify the configuration file name via a command line argument.')
    config_fname = argv[1]
#    print(sys.argv[2])
    # <2: เพื่อป้องกันข้อผิดพลาดจากการเข้าถึง sys.argv[1] เมื่อไม่มีการส่งอาร์กิวเมนต์เข้ามา
    # >2: เพื่อป้องกันข้อผิดพลาดจากการเข้าถึง sys.argv[2] เมื่อไม่มีการส่งอาร์กิวเมนต์ที่สองเข้ามา และตรวจสอบว่าค่านั้นตรงกับที่ต้องการหรือไม่


    # Initialize calculator, which diagonalizes mass/KE matrices, etc
    my_calculator = hubble_calculator(xin=config_fname)

    # Solve ODEs from tin to tfi
    my_calculator.solver()

#    ui,ux = my_calculator.solver()
#    plt.plot(ui,ux)

    # Save some information
#    my_calculator.output()

    # Make a plot
#    my_calculator.printout()


if __name__ == "__main__":
    main()
