from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle
import cycler

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 11
plt.rcParams["text.usetex"] = True
plt.rcParams["text.latex.preamble"] = r"""
\usepackage{physics}
\usepackage{siunitx}
%%% astronomical units
\DeclareSIUnit \au {au}
\DeclareSIUnit \jansky {Jk}
\DeclareSIUnit \ly {ly}
\DeclareSIUnit \parsec {pc}
\DeclareSIUnit \mag {mag}
\DeclareSIUnit \solarmass {\ensuremath{\mathit{M_{\odot}}}}
\DeclareSIUnit \yr {yr}
%%% cgs units
\DeclareSIUnit \erg {erg}
\DeclareSIUnit \cm {cm}
\DeclareSIUnit \gram {g}
"""

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

    # cosmos2020_cat = Table.read("COSMOS2020_photometric_only.fits", format="fits")
    # mask_bright = np.nonzero(
    #     np.logical_and.reduce(
    #         [
    #             (cosmos2020_cat["CFHT_u_FLUX_AUTO"] > 0.0),
    #             (cosmos2020_cat["HSC_g_FLUX_AUTO"] > 0.0),
    #             (cosmos2020_cat["HSC_r_FLUX_AUTO"] > 0.0),
    #             (cosmos2020_cat["HSC_i_FLUX_AUTO"] > 0.0),
    #             (cosmos2020_cat["HSC_z_FLUX_AUTO"] > 0.0),
    #             (cosmos2020_cat["HSC_y_FLUX_AUTO"] > 0.0),
    #             (cosmos2020_cat["UVISTA_J_FLUX_AUTO"] > 0.0),
    #             (cosmos2020_cat["UVISTA_H_FLUX_AUTO"] > 0.0),
    #             (cosmos2020_cat["UVISTA_Ks_FLUX_AUTO"] > 0.0),
    #         ]
    #     )
    # )
    # bright_gal_cat = cosmos2020_cat[mask_bright]
    # wavelengths = [3823.29, 4816.12, 6234.11, 7740.58, 9125.20, 9779.93, 12535.65, 16453.41, 21539.88]


# def get_column(name: str, col_id: int) -> float:
#     return bright_gal_cat[name][bright_gal_cat["ID"] == col_id][0]


# def get_fluxes(col_id: int) -> list[float]:
#     return np.ma.array(
#         [
#             get_column("CFHT_u_FLUX_AUTO", col_id),
#             get_column("HSC_g_FLUX_AUTO", col_id),
#             get_column("HSC_r_FLUX_AUTO", col_id),
#             get_column("HSC_i_FLUX_AUTO", col_id),
#             get_column("HSC_z_FLUX_AUTO", col_id),
#             get_column("HSC_y_FLUX_AUTO", col_id),
#             get_column("UVISTA_J_FLUX_AUTO", col_id),
#             get_column("UVISTA_H_FLUX_AUTO", col_id),
#             get_column("UVISTA_Ks_FLUX_AUTO", col_id),
#         ]
#     )


# def get_errors(col_id: int) -> list[float]:
#     return np.ma.array(
#         [
#             get_column("CFHT_u_FLUXERR_AUTO", col_id),
#             get_column("HSC_g_FLUXERR_AUTO", col_id),
#             get_column("HSC_r_FLUXERR_AUTO", col_id),
#             get_column("HSC_i_FLUXERR_AUTO", col_id),
#             get_column("HSC_z_FLUXERR_AUTO", col_id),
#             get_column("HSC_y_FLUXERR_AUTO", col_id),
#             get_column("UVISTA_J_FLUXERR_AUTO", col_id),
#             get_column("UVISTA_H_FLUXERR_AUTO", col_id),
#             get_column("UVISTA_Ks_FLUXERR_AUTO", col_id),
#         ]
#     )


def main():
    mode = darkmode(input("darkmode [d] or whitemode [Enter]?\n"))

    ### read data of 'COSMOS2020'
    cosmos2020_cat = Table.read("COSMOS2020_photometric_only.fits", format="fits")
    mask_bright = np.nonzero(
        np.logical_and.reduce(
            [
                (cosmos2020_cat["CFHT_u_FLUX_AUTO"] > 0.0),
                (cosmos2020_cat["HSC_g_FLUX_AUTO"] > 0.0),
                (cosmos2020_cat["HSC_r_FLUX_AUTO"] > 0.0),
                (cosmos2020_cat["HSC_i_FLUX_AUTO"] > 0.0),
                (cosmos2020_cat["HSC_z_FLUX_AUTO"] > 0.0),
                (cosmos2020_cat["HSC_y_FLUX_AUTO"] > 0.0),
                (cosmos2020_cat["UVISTA_J_FLUX_AUTO"] > 0.0),
                (cosmos2020_cat["UVISTA_H_FLUX_AUTO"] > 0.0),
                (cosmos2020_cat["UVISTA_Ks_FLUX_AUTO"] > 0.0),
            ]
        )
    )
    bright_gal_cat = cosmos2020_cat[mask_bright]


    def get_column(name: str, col_id: int) -> float:
        return bright_gal_cat[name][bright_gal_cat["ID"] == col_id][0]


    def get_fluxes(col_id: int) -> list[float]:
        return np.ma.array(
            [
                get_column("CFHT_u_FLUX_AUTO", col_id),
                get_column("HSC_g_FLUX_AUTO", col_id),
                get_column("HSC_r_FLUX_AUTO", col_id),
                get_column("HSC_i_FLUX_AUTO", col_id),
                get_column("HSC_z_FLUX_AUTO", col_id),
                get_column("HSC_y_FLUX_AUTO", col_id),
                get_column("UVISTA_J_FLUX_AUTO", col_id),
                get_column("UVISTA_H_FLUX_AUTO", col_id),
                get_column("UVISTA_Ks_FLUX_AUTO", col_id),
            ]
        )


    def get_errors(col_id: int) -> list[float]:
        return np.ma.array(
            [
                get_column("CFHT_u_FLUXERR_AUTO", col_id),
                get_column("HSC_g_FLUXERR_AUTO", col_id),
                get_column("HSC_r_FLUXERR_AUTO", col_id),
                get_column("HSC_i_FLUXERR_AUTO", col_id),
                get_column("HSC_z_FLUXERR_AUTO", col_id),
                get_column("HSC_y_FLUXERR_AUTO", col_id),
                get_column("UVISTA_J_FLUXERR_AUTO", col_id),
                get_column("UVISTA_H_FLUXERR_AUTO", col_id),
                get_column("UVISTA_Ks_FLUXERR_AUTO", col_id),
            ]
        )

    wavelengths = [3823.29, 4816.12, 6234.11, 7740.58, 9125.20, 9779.93, 12535.65, 16453.41, 21539.88]

    N = 6
    rng = np.random.default_rng()
    all_bright_gal_ids = [bright_gal["ID"] for bright_gal in bright_gal_cat]
    random_ids = rng.choice(all_bright_gal_ids, N, replace=False, shuffle=False)


    fluxes = {random_id: get_fluxes(random_id) for random_id in random_ids}
    errors = {random_id: get_errors(random_id) for random_id in random_ids}


    # Normalize
    for random_id in random_ids:
        scale_fac = 1.0 / np.nanmax(fluxes[random_id])
        fluxes[random_id] *= scale_fac
        errors[random_id] *= scale_fac


    # read data of filters
    DATA_DIR = "filters/CFHT_u.txt"
    CFHT_u_wavelength, CFHT_u_transmission = np.loadtxt(
        DATA_DIR,
        usecols=(0, 1),
        dtype=np.dtype([("CFHT_u_wavelength", float), ("CFHT_u_transmission", float)]),
        unpack=True,
    )

    DATA_DIR = "filters/HSC_g.txt"
    HSC_g_wavelength, HSC_g_transmission = np.loadtxt(
        DATA_DIR, usecols=(0, 1), dtype=np.dtype([("HSC_g_wavelength", float), ("HSC_g_transmission", float)]), unpack=True
    )

    DATA_DIR = "filters/HSC_r.txt"
    HSC_r_wavelength, HSC_r_transmission = np.loadtxt(
        DATA_DIR, usecols=(0, 1), dtype=np.dtype([("HSC_r_wavelength", float), ("HSC_r_transmission", float)]), unpack=True
    )

    DATA_DIR = "filters/HSC_i.txt"
    HSC_i_wavelength, HSC_i_transmission = np.loadtxt(
        DATA_DIR, usecols=(0, 1), dtype=np.dtype([("HSC_i_wavelength", float), ("HSC_i_transmission", float)]), unpack=True
    )

    DATA_DIR = "filters/HSC_y.txt"
    HSC_y_wavelength, HSC_y_transmission = np.loadtxt(
        DATA_DIR, usecols=(0, 1), dtype=np.dtype([("HSC_y_wavelength", float), ("HSC_y_transmission", float)]), unpack=True
    )

    DATA_DIR = "filters/HSC_z.txt"
    HSC_z_wavelength, HSC_z_transmission = np.loadtxt(
        DATA_DIR, usecols=(0, 1), dtype=np.dtype([("HSC_z_wavelength", float), ("HSC_z_transmission", float)]), unpack=True
    )

    DATA_DIR = "filters/UVISTA_H.txt"
    UVISTA_H_wavelength, UVISTA_H_transmission = np.loadtxt(
        DATA_DIR,
        usecols=(0, 1),
        dtype=np.dtype([("UVISTA_H_wavelength", float), ("UVISTA_H_transmission", float)]),
        unpack=True,
    )
    DATA_DIR = "filters/UVISTA_J.txt"
    UVISTA_J_wavelength, UVISTA_J_transmission = np.loadtxt(
        DATA_DIR,
        usecols=(0, 1),
        dtype=np.dtype([("UVISTA_J_wavelength", float), ("UVISTA_J_transmission", float)]),
        unpack=True,
    )

    DATA_DIR = "filters/UVISTA_Ks.txt"
    UVISTA_Ks_wavelength, UVISTA_Ks_transmission = np.loadtxt(
        DATA_DIR,
        usecols=(0, 1),
        dtype=np.dtype([("UVISTA_Ks_wavelength", float), ("UVISTA_Ks_transmission", float)]),
        unpack=True,
    )


    fig, ax = plt.subplots()


    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_cycle = cycle(colors)

    for random_id in random_ids:
        color = next(color_cycle)
        ax.plot(wavelengths, fluxes[random_id], color=color, label=f"ID={random_id}")
        ax.errorbar(wavelengths, fluxes[random_id], yerr=errors[random_id], fmt=".", capsize=5, color=color)

    # plt.plot(wavelengths, fluxes[random_id], label=f"ID={random_id}")
        # plt.errorbar(wavelengths, fluxes[random_id], yerr=errors[random_id], fmt=".", capsize=5)

    # plot filters
    plt.plot(CFHT_u_wavelength, CFHT_u_transmission, color="#6a006f", alpha=0.5, linewidth=1.0, linestyle="dotted")
    ax.fill_between(CFHT_u_wavelength, CFHT_u_transmission, color="#6a006f", alpha=0.15)

    plt.plot(HSC_g_wavelength, HSC_g_transmission, color="#4246C9", alpha=0.5, linewidth=1.0, linestyle="dotted")
    ax.fill_between(HSC_g_wavelength, HSC_g_transmission, color="#4246C9", alpha=0.15)

    plt.plot(HSC_r_wavelength, HSC_r_transmission, color="#3872C2", alpha=0.5, linewidth=1.0, linestyle="dotted")
    ax.fill_between(HSC_r_wavelength, HSC_r_transmission, color="#3872C2", alpha=0.15)

    plt.plot(HSC_i_wavelength, HSC_i_transmission, color="#2ABDAC", alpha=0.5, linewidth=1.0, linestyle="dotted")
    ax.fill_between(HSC_i_wavelength, HSC_i_transmission, color="#2ABDAC", alpha=0.15)

    plt.plot(HSC_y_wavelength, HSC_y_transmission, color="#2DCD69", alpha=0.5, linewidth=1.0, linestyle="dotted")
    ax.fill_between(HSC_y_wavelength, HSC_y_transmission, color="#2DCD69", alpha=0.15)

    plt.plot(HSC_z_wavelength, HSC_z_transmission, color="#FFBB30", alpha=0.5, linewidth=1.0, linestyle="dotted")
    ax.fill_between(HSC_z_wavelength, HSC_z_transmission, color="#FFBB30", alpha=0.15)

    plt.plot(UVISTA_H_wavelength, UVISTA_H_transmission, color="#FF9A30", alpha=0.5, linewidth=1.0, linestyle="dotted")
    ax.fill_between(UVISTA_H_wavelength, UVISTA_H_transmission, color="#FF9A30", alpha=0.15)

    plt.plot(UVISTA_J_wavelength, UVISTA_J_transmission, color="#FF7130", alpha=0.5, linewidth=1.0, linestyle="dotted")
    ax.fill_between(UVISTA_J_wavelength, UVISTA_J_transmission, color="#FF7130", alpha=0.15)

    plt.plot(UVISTA_Ks_wavelength, UVISTA_Ks_transmission, color="#E15A1F", alpha=0.5, linewidth=1.0, linestyle="dotted")
    ax.fill_between(UVISTA_Ks_wavelength, UVISTA_Ks_transmission, color="#E15A1F", alpha=0.15)

    plt.ylabel(r'Flux $[\si{\jansky}]$')
    plt.xlabel(r'$\lambda \, [\si{\angstrom}]$')
    plt.xlim(min(CFHT_u_wavelength), max(UVISTA_Ks_wavelength))
    plt.ylim(0.0)
    plt.legend(loc='upper left')
    plt.grid(True, linestyle='dashed')
    # plt.show()

    filename = 'wavelength-vs-flux'

    if mode==True:
        filename = filename + '_DARKMODE'

    # fig.savefig(filename + '.eps', format='eps', bbox_inches='tight')
    # fig.savefig(filename + '.pdf', format='pdf', bbox_inches='tight')
    fig.savefig(filename + '.png', format='png', bbox_inches='tight', dpi=300)
    # fig.savefig(filename + '.svg', format='svg', bbox_inches='tight')

if __name__ == "__main__":
    main()
