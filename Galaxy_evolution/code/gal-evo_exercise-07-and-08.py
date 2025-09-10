from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy import optimize as opt

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size']   = 11
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'''
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
'''


def binning(
        redshifts: npt.NDArray,
        bin_edges: npt.NDArray,
        samples_per_bin: int=10
        ) -> npt.NDArray:
    rng = np.random.default_rng()
    num_bins = len(bin_edges) - 1
    bin_indices = np.digitize(redshifts, bin_edges, right=False)
    return np.array([
        rng.choice(np.flatnonzero(bin_indices ==  i + 1), size=samples_per_bin) for i in range(num_bins)
        ])


def log10SFR_time_Eq10_fit(log_M_star, t, a0, a1, b0, b1, b2):
    return (a1 * t + b1) * log_M_star + b2 * log_M_star * log_M_star + (b0 + a0 * t)





def main():
    # Import data
    all_galaxies_table = Table.read('results_no_nebular_new_true.txt', format='ascii')

    log_ssfr = np.log10(all_galaxies_table['bayes.sfh.sfr10Myrs']/all_galaxies_table['bayes.stellar.m_star'])

    mask_quiescent = (log_ssfr < -10.5)
    mask_starforming = (log_ssfr > -10.5)

    quiescent_table = all_galaxies_table[mask_quiescent]
    starforming_table = all_galaxies_table[mask_starforming]

    # quiescent_table = Table.read('results_quiescent.txt', format='ascii')
    # starforming_table = Table.read('results_starforming.txt', format='ascii')

    # Plot age t [Gyr] against log10(SFR) [M_{*}/yr], divide 'best.universe.age' by 1000.0 to convert Myr to Gyr

    fig, ax = plt.subplots()

    plt.scatter(quiescent_table['best.universe.age']/1000.0, np.log10(quiescent_table['bayes.sfh.sfr10Myrs']),  color='tab:orange', marker='.', s=1.0, label='quiescent')
    plt.scatter(starforming_table['best.universe.age']/1000.0, np.log10(starforming_table['bayes.sfh.sfr10Myrs']), color='tab:blue', marker='.', s=1.0, label='starforming')
    plt.xlabel(r'age $t$ $[\si{\giga \yr}]$')
    plt.ylabel(r'$\log_{10}(SFR) \ [\si{\solarmass} / \si{\yr}$]')
    plt.legend(loc='upper left')
    plt.grid(True)
    # plt.show()

    figname = 'age-vs-log10SFR'
    # fig.savefig(figname + '.eps', format='eps', bbox_inches='tight')
    # fig.savefig(figname + '.pdf', format='pdf', bbox_inches='tight')
    fig.savefig(figname + '.png', format='png', bbox_inches='tight', dpi=250)
    # fig.savefig(figname + '.svg', format='svg', bbox_inches='tight')

    # Plot redshift z against log10(SFR) [M_{*}/yr]

    fig, ax = plt.subplots()

    plt.scatter(quiescent_table['best.universe.redshift'], np.log10(quiescent_table['bayes.sfh.sfr10Myrs']), color='tab:orange', marker='.', s=1.0, label='quiescent')
    plt.scatter(starforming_table['best.universe.redshift'], np.log10(starforming_table['bayes.sfh.sfr10Myrs']), color='tab:blue', marker='.', s=1.0, label='starforming')
    plt.xlabel(r'redshift $z$')
    plt.ylabel(r'$\log_{10}(SFR)$ [$\si{\solarmass} / \si{\yr}$]')
    plt.legend(loc='upper left')
    plt.grid(True)

    # plt.show()

    figname = 'redshift-vs-log10SFR'
    # fig.savefig(figname + '.eps', format='eps', bbox_inches='tight')
    # fig.savefig(figname + '.pdf', format='pdf', bbox_inches='tight')
    fig.savefig(figname + '.png', format='png', bbox_inches='tight', dpi=250)
    # fig.savefig(figname + '.svg', format='svg', bbox_inches='tight')


    # binning by redshift
    bin_edges = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8])
    starforming_redshift_binning = binning(starforming_table['best.universe.redshift'], bin_edges, 1000)
    all_galaxies_redshift_binning = binning(all_galaxies_table['best.universe.redshift'], bin_edges, 1000)

    # best-fit parameters according to Popesso et al. 2023, Table 1: https://arxiv.org/html/2403.06575v1#S4 or https://academic.oup.com/mnras/article/519/1/1526/6815739

    # all galaxies
    a1_guess_all_galaxies = 3.69
    a2_guess_all_galaxies = -3.81
    a3_guess_all_galaxies = 0.47
    b1_guess_all_galaxies = 11.91
    b2_guess_all_galaxies = -2.48
    b3_guess_all_galaxies = 0.44

    # starforming galaxies -- redshift fit
    # a1_guess_starforming = 3.47
    # a2_guess_starforming = -3.13
    # a3_guess_starforming = 0.56
    # b1_guess_starforming = 11.69
    # b2_guess_starforming = -1.66
    # b3_guess_starforming = 0.53

    # starforming galaxies -- Equation 10 time fit
    a0_guess_Eq10 = 0.20
    a1_guess_Eq10 = -0.034
    b0_guess_Eq10 = -26.134
    b1_guess_Eq10 = 4.722
    b2_guess_Eq10 = -0.1925

    # Equation 14 -- time fit
    a0_guess_Eq14 = 2.693
    a1_guess_Eq14 = -0.186
    a2_guess_Eq14 = 10.85
    a3_guess_Eq14 = -0.0729
    a4_guess_Eq14 = 0.99

    fig, ax = plt.subplots()

    sfr_SFG = starforming_table['bayes.sfh.sfr10Myrs']
    log10_sfr_SFG = np.log10(starforming_table['bayes.sfh.sfr10Myrs'])
    binned_log10_sfr_SFG = log10_sfr_SFG[starforming_redshift_binning]

    m_star_SFG = starforming_table['bayes.stellar.m_star']
    binned_m_star_SFG = m_star_SFG[starforming_redshift_binning]

    log10_m_star_SFG = np.log10(starforming_table['bayes.stellar.m_star'])
    binned_log10_m_star_SFG = np.log10(binned_m_star_SFG)

    redshift_SFG = starforming_table['best.universe.redshift']
    binned_redshift_SFG = redshift_SFG[starforming_redshift_binning]

    time_SFG = starforming_table['best.universe.age']
    binned_time_SFG = time_SFG[starforming_redshift_binning]

    cont_binned_m_star_SFG = np.linspace(np.min(binned_m_star_SFG), np.max(binned_m_star_SFG), 1000)
    log10_cont_binned_m_star_SFG = np.log10(cont_binned_m_star_SFG)

    cont_m_star_SFG = np.logspace(np.log10(np.min(m_star_SFG)), np.log10(np.max(m_star_SFG)), 1000)
    log10_cont_m_star_SFG = np.log10(cont_m_star_SFG)

    c = binned_redshift_SFG
    c = np.repeat(bin_edges[:-1].reshape(-1, 1), len(starforming_redshift_binning[0]), axis=1)

    # end = 0.4 * len(bin_edges)
    # offset = np.repeat(np.arange(0, end - 0.5, 0.4).reshape(-1, 1), len(starforming_redshift_binning[0]), axis=1)

    # plt.scatter(log10_m_star_SFG[starforming_redshift_binning], log10_sfr_SFG[starforming_redshift_binning] + offset, marker='.', s=14.0, c=c, cmap="hsv")
    plt.scatter(log10_m_star_SFG[starforming_redshift_binning], log10_sfr_SFG[starforming_redshift_binning], marker='.', s=14.0, c=c, cmap="hsv")
    # plt.show()

    times = np.linspace(1.0, 13.0, 13)
    # [print(zi) for zi in redshift_bin_mids]

    # ax.set_prop_cycle('color', plt.cm.hsv(np.linspace(0, 1, len(bin_edges))))
    ax.set_prop_cycle('color', plt.cm.hsv(np.linspace(0, 1, len(times))))

    # [plt.plot(log10_cont_binned_m_star_SFG, log10SFR_redshift_fit(cont_binned_m_star_SFG, zi, a1_guess_starforming, a2_guess_starforming, a3_guess_starforming, b1_guess_starforming, b2_guess_starforming, b3_guess_starforming), lw=1.0) for zi in redshift_bin_mids]
    [plt.plot(log10_cont_binned_m_star_SFG, log10SFR_time_Eq10_fit(log10_cont_binned_m_star_SFG, ti, a0_guess_Eq10, a1_guess_Eq10, b0_guess_Eq10, b1_guess_Eq10, b2_guess_Eq10), lw=1.0, linestyle='dashed') for ti in times]
    [plt.plot(log10_cont_binned_m_star_SFG, log10SFR_time_Eq14_fit(cont_binned_m_star_SFG, ti, a0_guess_Eq14, a1_guess_Eq14, a2_guess_Eq14, a3_guess_Eq14, a4_guess_Eq14), lw=1.0) for ti in times]
#        plt.plot(log_cont_binned_m_star, log10sfr_fit, cmap="hsv")
    # colors = [matplotlib.cm.hsv(x) for x in redshift_bin_mids]
    # print("colors = ", colors)

    # for i in range(len(redshift_bin_mids)):
    #     z = redshift_bin_mids[i]
    #     log10sfr_fit = log10SFR_redshift_fit(cont_binned_m_star, z, a1_guess_starforming, a2_guess_starforming, a3_guess_starforming, b1_guess_starforming, b2_guess_starforming, b3_guess_starforming)
# #        plt.plot(log_cont_binned_m_star, log10sfr_fit, cmap="hsv")
    #     plt.plot(log_cont_binned_m_star, log10sfr_fit, lw=1.0)

    plt.xlabel(r'$\log_{10}(M_{*}) \ [\si{\solarmass}]$ ')
    plt.ylabel(r'$\log_{10}(SFR) \ [\si{\solarmass} / \si{\yr}]$ ')
    plt.title(r'Starforming Galaxies')
    # plt.xlim(8.6, 12.0)
    # plt.ylim(-1.0, 2.9)
    plt.grid(True)

    figname = 'TIME_log10_m_star_vs_log10_sfr_binned_SFG.png'
    # fig.savefig(figname + '.eps', format='eps', bbox_inches='tight')
    # fig.savefig(figname + '.pdf', format='pdf', bbox_inches='tight')
    fig.savefig(figname + '.png', format='png', bbox_inches='tight', dpi=250)
    # fig.savefig(figname + '.svg', format='svg', bbox_inches='tight')

    fig, ax = plt.subplots()

    sfr_ALL = all_galaxies_table['bayes.sfh.sfr10Myrs']
    log10_sfr_ALL = np.log10(all_galaxies_table['bayes.sfh.sfr10Myrs'])
    binned_log10_sfr_ALL = log10_sfr_ALL[all_galaxies_redshift_binning]

    m_star_ALL = all_galaxies_table['bayes.stellar.m_star']
    binned_m_star_ALL = m_star_ALL[all_galaxies_redshift_binning]

    log10_m_star_ALL = np.log10(all_galaxies_table['bayes.stellar.m_star'])
    binned_log10_m_star_ALL = np.log10(binned_m_star_ALL)

    redshifts_ALL = all_galaxies_table['best.universe.redshift']
    binned_redshift_ALL = redshifts_ALL[all_galaxies_redshift_binning]

    time_ALL = all_galaxies_table['best.universe.age']
    binned_time_ALL = time_ALL[all_galaxies_redshift_binning]

    cont_m_star_ALL = np.logspace(np.log10(np.min(m_star_ALL)), np.log10(np.max(m_star_ALL)), 1000)
    log10_cont_m_star_ALL = np.log10(cont_m_star_ALL)

    c = binned_redshift_ALL
    # c = np.repeat(bin_edges[:-1].reshape(-1, 1), len(all_galaxies_redshift_binning[0]), axis=1)

    # ax.set_prop_cycle('color', plt.cm.hsv(np.linspace(0, 1, len(bin_edges))))

    # plt.scatter(log10_m_star_ALL[all_galaxies_redshift_binning], log10_sfr_ALL[all_galaxies_redshift_binning], marker='.', s=14.0, c=c, cmap="hsv")
    # [plt.plot(log10_cont_m_star_ALL, log10SFR_time_fit(cont_m_star_ALL, ti, a1_guess_all_galaxies, a2_guess_all_galaxies, a3_guess_all_galaxies, b1_guess_all_galaxies, b2_guess_all_galaxies, b3_guess_all_galaxies), lw=1.0) for ti in times]
    # [plt.plot(log10_cont_m_star_ALL, log10SFR_time_fit(cont_m_star_ALL, ti, a0_guess_all_galaxies, a1_guess_all_galaxies, b0_guess_all_galaxies, b1_guess_all_galaxies, b2_guess_all_galaxies), lw=1.0) for ti in times]

    # plt.xlabel(r'$\log_{10}(M_{*}) \ [M_{\odot}]$ ')
    # plt.ylabel(r'$\log_{10}(SFR) \ [M_{\odot} / \SI{}{\yr}]$ ')
    # plt.title(r'All galaxies')
    # plt.xlim(9.0, 11.5)
    # plt.ylim(-0.1, 2.9)
    # plt.grid(True)


    # fig.savefig('TIME_log10_m_star_-vs_log10_sfr_binned_ALL.png', format='png', bbox_inches='tight', dpi=250)

    # for i in range(len(bin_edges) - 1):
    #     print("starforming_table['best.universe.redshift'][starforming_redshift_binning][i] = ", starforming_table['best.universe.redshift'][starforming_redshift_binning][i])
    #     print("starforming_table['bayes.stellar.m_star'][starforming_redshift_binning][i] = ", starforming_table['bayes.stellar.m_star'][starforming_redshift_binning][i])
    #     SFR = starforming_table['bayes.sfh.sfr10Myrs'][starforming_redshift_binning][i]
    #     z = starforming_table['best.universe.redshift'][starforming_redshift_binning][i]
    #     M_star = starforming_table['bayes.stellar.m_star'][starforming_redshift_binning][i]

    #     log10SFR_fit = log10SFR_redshift_fit(M_star, z, a1_guess_starforming, a2_guess_starforming, a3_guess_starforming, b1_guess_starforming, b2_guess_starforming, b3_guess_starforming)
    #     print("log10SFR_redshift_fit = ", log10SFR_redshift_fit(M_star, z, a1_guess_starforming, a2_guess_starforming, a3_guess_starforming, b1_guess_starforming, b2_guess_starforming, b3_guess_starforming))
    #     print("np.log10(SFR) = ", np.log10(SFR))
    #     log_m_star = np.log10(starforming_table['bayes.stellar.m_star'])
    #     plt.scatter(log_m_star[starforming_redshift_binning], np.log10(SFR))





# ### FITTING ###
# ### ------- ###

# Starforming Galaxies
# --------------------
    fig, ax = plt.subplots()

    # Global fit..
    p0_guess = [a0_guess_Eq10,
                a1_guess_Eq10,
                b0_guess_Eq10,
                b1_guess_Eq10,
                b2_guess_Eq10]

    mass_and_time = np.array([log10_m_star_SFG, time_SFG])
    popt, pcov = opt.curve_fit(lambda x, *args: log10SFR_time_Eq10_fit(*x, *args), mass_and_time, log10_sfr_SFG, p0=p0_guess, maxfev=50000)

    print("-- estimated parameters --")
    print("*popt = ", *popt)
    print("a0_starforming = ", popt[0])
    print("a1_starforming = ", popt[1])
    print("b0_starforming = ", popt[2])
    print("b1_starforming = ", popt[3])
    print("b2_starforming = ", popt[4])

    ax.set_prop_cycle('color', plt.cm.hsv(np.linspace(0, 1, len(bin_edges))))
    # redshift_bin_mids = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.25, 4.25, 5.25])
    # times = np.linspace(1.0, 13.0, 13)
    for i in range(len(bin_edges) - 1):

        ti = binned_time_SFG[i]

        fitted_binned_log10_sfr = log10SFR_time_Eq10_fit(log10_cont_m_star_SFG, np.mean(binned_time_SFG[i]), *popt)
        # fitted_binned_log10_sfr = log10SFR_time_Eq10_fit(np.log10(cont_m_star), np.mean(binned_time_SFG[i]), *p0_guess)

        # plt.plot(log10_cont_m_star_SFG, fitted_binned_log10_sfr + offset[i], linewidth=1.0)
        plt.plot(log10_cont_m_star_SFG, fitted_binned_log10_sfr, linewidth=1.0)
        plt.xlabel(r'$\log_{10}(M_{*}) \ [\si{\solarmass}]$ ')
        plt.ylabel(r'$\log_{10}(SFR) \ [\si{\solarmass} / \si{\yr}]$')
        plt.title(r'My Fit of Starforming Galaxies')
        plt.grid(True)


    # plt.scatter(log10_m_star_SFG[starforming_redshift_binning], log10_sfr_SFG[starforming_redshift_binning] + offset, marker='.', s=14.0, c=c, cmap="hsv")
    plt.scatter(log10_m_star_SFG[starforming_redshift_binning], log10_sfr_SFG[starforming_redshift_binning], marker='.', s=14.0, c=c, cmap="hsv")
    # plt.xlim(9.0, 12.0)
    # plt.ylim(-1.0, 3.0)
    # fig.savefig('SFG_binned_log10_m_star-vs-fitted_binned_log10_sfr.png', format='png', bbox_inches='tight', dpi=250)
    fig.savefig('SFG_binned_m_star-vs-fitted_binned_log10_sfr.png', format='png', bbox_inches='tight', dpi=250)


# All Galaxies
# ------------
    fig, ax = plt.subplots()

    # Global fit..
    p0_guess = [a0_guess_Eq10,
                a1_guess_Eq10,
                b0_guess_Eq10,
                b1_guess_Eq10,
                b2_guess_Eq10]

    mass_and_time = np.array([log10_m_star_ALL, time_ALL])
    popt, pcov = opt.curve_fit(lambda x, *args: log10SFR_time_Eq10_fit(*x, *args),  mass_and_time, log10_sfr_ALL, p0=p0_guess, maxfev=50000)

    print("-- estimated parameters --")
    print("*popt = ", *popt)
    print("a0_Eq10 = ", popt[0])
    print("a1_Eq10 = ", popt[1])
    print("b0_Eq10 = ", popt[2])
    print("b1_Eq10 = ", popt[3])
    print("b2_Eq10 = ", popt[4])

    ax.set_prop_cycle('color', plt.cm.hsv(np.linspace(0, 1, len(bin_edges))))
    # redshift_bin_mids = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.25, 4.25, 5.25])
    # times = np.linspace(1.0, 13.0, 13)
    for i in range(len(bin_edges) - 1):

        ti = binned_time_ALL[i]

        fitted_binned_log10_sfr = log10SFR_time_Eq10_fit(log10_cont_m_star_ALL, np.mean(binned_time_ALL[i]), *popt)

        # plt.plot(log10_cont_m_star_ALL, fitted_binned_log10_sfr + offset[i], linewidth=1.0)
        plt.plot(log10_cont_m_star_ALL, fitted_binned_log10_sfr, linewidth=1.0)
        plt.xlabel(r'$\log_{10}(M_{*}) \ [\si{\solarmass}]$ ')
        plt.ylabel(r'$\log_{10}(SFR) \ [\si{\solarmass} / \si{\yr}]$')
        plt.title(r'My Fit for All Galaxies')
        # plt.xlim(5.0, 17.5)
        # plt.ylim(0.0, 12.0)
        plt.grid(True)


    plt.scatter(log10_m_star_ALL[all_galaxies_redshift_binning], log10_sfr_ALL[all_galaxies_redshift_binning], marker='.', s=14.0, c=c, cmap="hsv")
    # plt.scatter(log10_m_star_ALL[all_galaxies_redshift_binning], log10_sfr_ALL[all_galaxies_redshift_binning] + offset, marker='.', s=14.0, c=c, cmap="hsv")
    # plt.xlim(9.0, 12.0)
    # plt.ylim(-1.0, 3.0)
    # fig.savefig('ALL_binned_log10_m_star-vs-fitted_binned_log10_sfr.png', format='png', bbox_inches='tight', dpi=250)
    fig.savefig('ALL_binned_log10_m_star-vs-fitted_binned_log10_sfr.png', format='png', bbox_inches='tight', dpi=250)


if __name__ == "__main__":
    main()
