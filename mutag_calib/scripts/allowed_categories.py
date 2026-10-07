""" Define categories, that are to be processed by the following scripts:

- ALLOWED_CATEGORIES
    - run_fit_results.py
    - run_all_combine_plots.py
This is just a list of categories that are to be processed.
Any item in the list that is not in the datacards is just ignored. Therefore, it doesn't hurt to have too many categories active.


- ALLOWED_CATEGORIES_SF_PLOT
    - make_SFs_plots.py
This is a dictionary of different binnings for SF combinations. Every dictionary entry represents one set of SFs that are to be combined into one set.
These should ideally be mutually exclusive. No checks are performed in that regard.
"""

ALLOWED_CATEGORIES = {
    # HHbbtt / HHbbgg categories (from upstream)
    "msd-80to170_Pt-300to350_particleNet_XbbVsQCD-HHbbtt",
    "msd-80to170_Pt-350to425_particleNet_XbbVsQCD-HHbbtt",
    "msd-80to170_Pt-425toInf_particleNet_XbbVsQCD-HHbbtt",
    "msd-30toInf_Pt-300to350_particleNet_XbbVsQCD-HHbbgg",
    "msd-30toInf_Pt-350to425_particleNet_XbbVsQCD-HHbbgg",
    "msd-30toInf_Pt-425toInf_particleNet_XbbVsQCD-HHbbgg",
    "msd-30toInf_Pt-300to350_globalParT3_XbbVsQCD-HHbbgg",
    "msd-30toInf_Pt-350to425_globalParT3_XbbVsQCD-HHbbgg",
    "msd-30toInf_Pt-425toInf_globalParT3_XbbVsQCD-HHbbgg",
    "msd-30toInf_Pt-300to400_globalParT3_XbbVsQCD-HHbbgg",
    "msd-30toInf_Pt-400to450_globalParT3_XbbVsQCD-HHbbgg",
    "msd-30toInf_Pt-450toInf_globalParT3_XbbVsQCD-HHbbgg",

    # HHbbbb categories
    # "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb",
    # "msd-100to150_Pt-250toInf_particleNet_XbbVsQCD-HHbbbb",
    # "msd-100to150_Pt-250to300_particleNet_XbbVsQCD-HHbbbb",
    # "msd-100to150_Pt-300to350_particleNet_XbbVsQCD-HHbbbb",
    # "msd-100to150_Pt-350to425_particleNet_XbbVsQCD-HHbbbb",
    # "msd-100to150_Pt-425toInf_particleNet_XbbVsQCD-HHbbbb",
    # "msd-100to150_Pt-250to300_globalParT3_XbbVsQCD-HHbbbb",
    # "msd-100to150_Pt-300to350_globalParT3_XbbVsQCD-HHbbbb",
    # "msd-100to150_Pt-350to425_globalParT3_XbbVsQCD-HHbbbb",
    # "msd-100to150_Pt-425toInf_globalParT3_XbbVsQCD-HHbbbb",

    # "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_WP_MP",
    # "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_WP_MP",
    # "msd-50to200_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_WP_MP",
    # "msd-50to200_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_WP_MP",
    # "msd-100to150_Pt-250to330_globalParT3_XbbVsQCD-HHbbbb_LP",
    # "msd-100to150_Pt-330to420_globalParT3_XbbVsQCD-HHbbbb_LP",
    # "msd-100to150_Pt-420toInf_globalParT3_XbbVsQCD-HHbbbb_LP",
    # "msd-100to150_Pt-250to330_globalParT3_XbbVsQCD-HHbbbb_MP",
    # "msd-100to150_Pt-330to420_globalParT3_XbbVsQCD-HHbbbb_MP",
    # "msd-100to150_Pt-420toInf_globalParT3_XbbVsQCD-HHbbbb_MP",
    # "msd-100to150_Pt-250to330_globalParT3_XbbVsQCD-HHbbbb_HP",
    # "msd-100to150_Pt-330to420_globalParT3_XbbVsQCD-HHbbbb_HP",
    # "msd-100to150_Pt-420toInf_globalParT3_XbbVsQCD-HHbbbb_HP",
    # "msd-100to150_Pt-250to330_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    # "msd-100to150_Pt-330to420_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    # "msd-100to150_Pt-420toInf_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",

    # "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_LP",
    # "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_LP",
    # "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_MP",
    # "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_MP",
    # "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_HP",
    # "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_HP",

    # # "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_WP_MP",
    # # "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_WP_MP",
    # # "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_WP_HP",
    # # "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_WP_HP",
    # "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    # "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",

    # "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_LP",
    # "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_MP",
    # "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_HP",
    # # "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_WP_MP",
    # # "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_WP_HP",
    # "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",

    "msd-100to150_Pt-250to300_globalParT3_XbbVsQCD-HHbbbb_LP",
    "msd-100to150_Pt-300to400_globalParT3_XbbVsQCD-HHbbbb_LP",
    "msd-100to150_Pt-400to500_globalParT3_XbbVsQCD-HHbbbb_LP",
    "msd-100to150_Pt-500toInf_globalParT3_XbbVsQCD-HHbbbb_LP",
    "msd-100to150_Pt-250to300_globalParT3_XbbVsQCD-HHbbbb_MP",
    "msd-100to150_Pt-300to400_globalParT3_XbbVsQCD-HHbbbb_MP",
    "msd-100to150_Pt-400to500_globalParT3_XbbVsQCD-HHbbbb_MP",
    "msd-100to150_Pt-500toInf_globalParT3_XbbVsQCD-HHbbbb_MP",
    "msd-100to150_Pt-250to300_globalParT3_XbbVsQCD-HHbbbb_HP",
    "msd-100to150_Pt-300to400_globalParT3_XbbVsQCD-HHbbbb_HP",
    "msd-100to150_Pt-400to500_globalParT3_XbbVsQCD-HHbbbb_HP",
    "msd-100to150_Pt-500toInf_globalParT3_XbbVsQCD-HHbbbb_HP",
    "msd-100to150_Pt-250to300_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    "msd-100to150_Pt-300to400_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    "msd-100to150_Pt-400to500_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    "msd-100to150_Pt-500toInf_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
}

ALLOWED_CATEGORIES_SF_PLOT = {
    # "single_bin": ["msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_WP_MP"],
    # "multi_purities": ["msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_LP",
    #                    "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_MP",
    #                    "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_HP",
    #                    "msd-100to150_Pt-250toInf_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    #                    ],
    # "multi_purities_pt_bins": [
    #                    "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_LP",
    #                    "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_LP",
    #                    "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_MP",
    #                    "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_MP",
    #                    "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_HP",
    #                    "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_HP",
    #                    "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    #                    "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    #                    ],
    "multi_purities_4pt_bins": [
                "msd-100to150_Pt-250to300_globalParT3_XbbVsQCD-HHbbbb_LP",
                "msd-100to150_Pt-300to400_globalParT3_XbbVsQCD-HHbbbb_LP",
                "msd-100to150_Pt-400to500_globalParT3_XbbVsQCD-HHbbbb_LP",
                "msd-100to150_Pt-500toInf_globalParT3_XbbVsQCD-HHbbbb_LP",
                "msd-100to150_Pt-250to300_globalParT3_XbbVsQCD-HHbbbb_MP",
                "msd-100to150_Pt-300to400_globalParT3_XbbVsQCD-HHbbbb_MP",
                "msd-100to150_Pt-400to500_globalParT3_XbbVsQCD-HHbbbb_MP",
                "msd-100to150_Pt-500toInf_globalParT3_XbbVsQCD-HHbbbb_MP",
                "msd-100to150_Pt-250to300_globalParT3_XbbVsQCD-HHbbbb_HP",
                "msd-100to150_Pt-300to400_globalParT3_XbbVsQCD-HHbbbb_HP",
                "msd-100to150_Pt-400to500_globalParT3_XbbVsQCD-HHbbbb_HP",
                "msd-100to150_Pt-500toInf_globalParT3_XbbVsQCD-HHbbbb_HP",
                "msd-100to150_Pt-250to300_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
                "msd-100to150_Pt-300to400_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
                "msd-100to150_Pt-400to500_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
                "msd-100to150_Pt-500toInf_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
        ],
    # "multi_purities_3pt_bins": [
    #                     "msd-100to150_Pt-250to330_globalParT3_XbbVsQCD-HHbbbb_LP",
    #                     "msd-100to150_Pt-330to420_globalParT3_XbbVsQCD-HHbbbb_LP",
    #                     "msd-100to150_Pt-420toInf_globalParT3_XbbVsQCD-HHbbbb_LP",
    #                     "msd-100to150_Pt-250to330_globalParT3_XbbVsQCD-HHbbbb_MP",
    #                     "msd-100to150_Pt-330to420_globalParT3_XbbVsQCD-HHbbbb_MP",
    #                     "msd-100to150_Pt-420toInf_globalParT3_XbbVsQCD-HHbbbb_MP",
    #                     "msd-100to150_Pt-250to330_globalParT3_XbbVsQCD-HHbbbb_HP",
    #                     "msd-100to150_Pt-330to420_globalParT3_XbbVsQCD-HHbbbb_HP",
    #                     "msd-100to150_Pt-420toInf_globalParT3_XbbVsQCD-HHbbbb_HP",
    #                     "msd-100to150_Pt-250to330_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    #                     "msd-100to150_Pt-330to420_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    #                     "msd-100to150_Pt-420toInf_globalParT3_XbbVsQCD-HHbbbb_WP_VHP",
    #                    ],
    # "multi_pt": [
    #     "msd-100to150_Pt-250to350_globalParT3_XbbVsQCD-HHbbbb_WP_MP",
    #     "msd-100to150_Pt-350toInf_globalParT3_XbbVsQCD-HHbbbb_WP_MP",
    # ],
}
