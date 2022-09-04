
	0. GENERAL INFORMATION

ADV measurements conducted in three cross sections:

- Downstream
- Middle
- Upstream

under 5 flow conditions:

ID		V1	V2	V3	V4	V5
q (m2/s)	0.09	0.13	0.19	0.23	0.26

All measurements were sampled every 2 cm. An additional measurement was conducted in the "Downstream" cross-section using a finer spacing (1 cm): "Downstream-fine"

	1. RAW DATA

The following code can be used to access the raw data: "ADV_profile_iterator.py"

The following ad-hoc libraries are also required: "faux_input", "faux_filter".

When running ADV_profile_iterator, the code runs over all data folders looking for a txt file named "info.txt", containing the following keywords (experiment dependent):

q_step	Zmin(m)	Dz(m)	Hmax(m)	ks(m)
3	0.000	0.02	0.278	6.7E-3

here explained:
- q_step: # in the flow ID (for instance, 3 for V3).
- Zmin(m): start of the flow sampling (always 0).
- Dz(m): vertical resolution in the flow sampling (0.02 m, expcept for "Downstream-fine").
- Hmax(m): flow depth (0.278 m for all experiments).
- ks(m): sand roughness (estimated at 6.7E-3 m).

If a "info.txt" is found inside a folder (for instance, at "Downstream/V3"), that profile is analysed. Relevant plots are produced including mean and turbulent quantities and txt files are generated saving the profile final results.

I recommend to only have one "info.txt" file in all folders (not tested with several).

	2. PROCESSED DATA

All analysed profiles are in this folder, which includes the following subfolders with processed ADV data:

- Downstream: data on the downstream cross-section with measuring spacing 2 cm.
- Downstream-fine: data on the downstream cross-section with measuring spacing 1 cm.
- Middle: data on the middle cross-section with measuring spacing 2 cm.
- Upstream: data on the upstream cross-section with measuring spacing 2 cm.

Besides, the following subfolders include digitsed data from literature: "_liter-data stream fluctuations", "_liter-data normal fluctuations", with streamwise and normalwise (normal to bed) fluctuations.

The code(s) [Plot_figure_Vx_uxux_uzuz] can be run to generate figures showing mean velocity, streamwise fluctuation and normalwise fluctuation in three cross-sections: downstream [_DS], middle [_MID], upstream [_US].
