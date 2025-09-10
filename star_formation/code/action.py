### astropy package -- https://docs.astropy.org/en/stable/ ###
from astropy import constants as const
from astropy import units as u
from astropy.units import imperial
### matplotlib package -- https://matplotlib.org/stable/index.html ###
import matplotlib.pyplot as plt
### numpy package -- https://numpy.org/doc/stable/ ###
import numpy as np
### scipy package -- https://docs.scipy.org/doc/scipy/ ###
# import scipy.constants as const
from scipy import optimize as opt                      #   for optimization and fit -- https://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size']   = 11
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'''
\usepackage{physics}
\usepackage{siunitx}

\AtBeginDocument{\RenewCommandCopy\qty\SI}

%%% Astronomical units
\DeclareSIUnit \au {au}
\DeclareSIUnit \ly {ly}
\DeclareSIUnit \parsec {pc}
\DeclareSIUnit \mag {mag}
\DeclareSIUnit \Jy {Jy}
% 1 Jy = 1e-26 W m^{-2} Hz^{-1} = 1e-23 erg s^{-1} cm^{-2} Hz^{-1}
\DeclareSIUnit \yr {yr}
%%% CGS units
\DeclareSIUnit \erg {erg}
\DeclareSIUnit \cm {cm}
'''


def darkmode(colormode: str):
    if colormode=="d":
        plt.rcParams["backend"] = "Qt5Agg"
        plt.rcParams["text.color"] = "#FFFFFF"
        plt.rcParams["axes.facecolor"] = "#131A24"
        plt.rcParams["axes.edgecolor"] = "#FFFFFF"
        plt.rcParams["axes.grid"] = True
        plt.rcParams["axes.labelcolor"] = "#FFFFFF"
        plt.rcParams["axes.formatter.use_mathtext"] = True
        plt.rcParams["axes.prop_cycle"] = cycler.cycler('color', ['#007FFF', '#00FF7F', '#FF007F', '#10BCEB', '#7F77FF', '#FFF777'])
        plt.rcParams["xtick.color"] = "#FFFFFF"
        plt.rcParams["ytick.color"] = "#FFFFFF"
        plt.rcParams["grid.color"] = "#b0b0b0"
        plt.rcParams["grid.linestyle"] = "dashed"
        plt.rcParams["figure.facecolor"] = "#131A24"
        plt.rcParams["figure.edgecolor"] = "#131A24"
        return True
    else:
        return False


def flux(F0, m):
    return F0 * np.power(10.0, -m/2.5)


def flux_err(F0, m, m_err):
    return np.abs((-1.0/2.5) * np.log(10.0) * F0 * np.power(10.0, -m/2.5) * m_err)


def fit(x, a, b):
    return a * x + b


def main():
    mode = darkmode(input("darkmode [d] or whitemode [Enter]?\n"))
    ### IMPORT DATA ###
    ### ----------- ###
    DATA_DIR = 'table.txt'
    ID, X_pos, Y_pos, Mag_1, Mag_err_1, Mag_2, Mag_err_2, Mag_3, Mag_err_3, Mag_4, Mag_err_4 = np.loadtxt(DATA_DIR,
        usecols=(0,1,2,3,4,5,6,7,8,9,10),
        dtype=np.dtype([('ID', int),
                        ('X_pos', float),
                        ('Y_pos', float),
                        ('Mag_1', float),
                        ('Mag_err_1', float),
                        ('Mag_2', float),
                        ('Mag_err_2', float),
                        ('Mag_3', float),
                        ('Mag_err_3', float),
                        ('Mag_4', float),
                        ('Mag_err_4', float)
                        ]),
        comments='#',
        unpack=True)

    print("IDs = ", ID)
    # Create a dictionary with the variables
    magnitudes = {f"Mag_{x}": locals()[f"Mag_{x}"] for x in range(1, 5)}
    magnitude_errors = {f"Mag_err_{x}": locals()[f"Mag_err_{x}"] for x in range(1, 5)}


    DATA_DIR = 'wavelength-flux.txt'
    Band, wavelength, flux_density = np.loadtxt(DATA_DIR,
        usecols=(0,1,2),
        dtype=np.dtype([('Band', str),
                        ('wavelength', float),
                        ('flux_density', float)
                        ]),
        comments='#',
        unpack=True)

    # some guesses for the slope of linear fit in [micro meter/Jansky]
    # a_parameter_guess = [
    #         (8.0 - 3.6)/(0.07 - 0.068),
    #         (8.0 - 3.6)/(0.06 - 0.250),
    #         (8.0 - 3.6)/(0.165 - 0.25),
    #         (8.0 - 3.6)/(0.05 - 0.13),
    #         (8.0 - 4.0)/(0.355 - 0.335),
    #         (8.0 - 4.0)/(0.36 - 0.28),
    #         (8.0 - 3.6)/(0.36 - 0.05),
    #         (8.0 - 3.6)/(0.04 - 0.029),
    #         (8.0 - 4.5)/(0.271 - 0.251),
    #         (8.0 - 3.6)/(0.24 - 0.17)
    #         ]

    print("=========")
    print("Band = ", Band)
    print("wavelength = ", wavelength)
    print("flux_density = ", flux_density)
    print("=========")

    for i in range(len(ID)):
        fig, ax = plt.subplots()

        print("\nID = ", ID[i])
        print("-----------")

        wl_LIST = []
        F_LIST = []
        m_LIST = []
        m_err_LIST = []

        for j in range(4):
            wl = wavelength[3 + j]
            F0 = flux_density[3 + j]
            m = magnitudes[f"Mag_{j+1}"][i]
            m_err = magnitude_errors[f"Mag_err_{j+1}"][i]

            flux_error = flux_err(F0, m, m_err)

            wl_LIST.append(wl)
            F_LIST.append(flux(F0,m))
            m_LIST.append(m)

            print("lambda = ", wl)
            print("F0 = ", F0)
            print(f"Mag_{j+1} = ", m)
            print("F = ", flux(F0,m))
            # print(f"Mag_err_{j+1} = ", magnitude_errors[f"Mag_err_{j+1}"][i])
            print(f"flux_err_{j+1} = ", flux_error)

            # plt.scatter(wl, flux(F0, m), s=4.0, label=f'IRAC{j + 1}')
            plt.errorbar(wl, flux(F0, m), yerr=flux_error, capsize=4.0, label=f'IRAC{j + 1}')
            # plt.errorbar(wl, flux(F0, m), yerr=flux_error, capsize=4.0)

        # a = a_parameter_guess[i]

        popt, pcov = opt.curve_fit(fit, wl_LIST, F_LIST, p0=[0.0,0.1])
        # print("popt = ", popt)
        print("a_popt = ", popt[0])
        print("a_eq = ", flux(F0,m)/(flux(F0,m) - wl))
        wl_cont = np.linspace(min(wl_LIST), max(wl_LIST), 1000)
        # plt.plot(wl, flux(F0,m), label=f'ID = {ID[i]}')
        plt.plot(wl, flux(F0,m))
        F_fit = fit(wl_cont, *popt)
        plt.plot(wl_cont, F_fit, label=f'linear fit: ${popt[0]:.4f} x + {popt[1]:.4f}$', color='purple')


        plt.xlabel(r'$\lambda$ in [$\si{\micro \meter}$]')
        plt.ylabel(r'$F$ in [$\si{\Jy}$]')
        plt.title(rf'ID = {ID[i]}')
        plt.legend(loc='upper right')
        plt.grid(True)

        fig.savefig(f'ID_{ID[i]}.png', format='png', bbox_inches='tight', dpi=250)
    # fig.savefig('STAR.png', format='png', bbox_inches='tight', dpi=250)




    ### PLOT ###
    ### ---- ###


    # fig, ax = plt.subplots()
    # plt.plot(X, Y, color='tab:blue', linewidth=1.0, label='X vs. Y')
    # plt.errorbar(X, Y, xerr=X_error, yerr=Y_error, fmt='.', ecolor='tab:red', elinewidth=1.0, label='error bars')
    # plt.xlabel(r'$X$ in \SI{}{\nano \meter}')
    # plt.ylabel(r'$Y$ in \SI{}{\nano \meter}')
    # plt.title(r'Some Title')
    # plt.legend(loc='upper right')
    # plt.grid(True)

    # fig.savefig('X-vs-Y.eps', format='eps', bbox_inches='tight')
    # fig.savefig('X-vs-Y.pdf', format='pdf', bbox_inches='tight')
    # fig.savefig('X-vs-Y.png', format='png', bbox_inches='tight', dpi=250)
    # fig.savefig('X-vs-Y.svg', format='svg', bbox_inches='tight')



if __name__ == "__main__":
    main()
