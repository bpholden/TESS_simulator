SHELL = /bin/sh
INSTALL = install -c
MKDIR = install -d 

destdir = ..

srcdir = .

fitdir = ../../SystPy/

SRC = ExposureCalculations.py \
	Generate_Errors.py \
	NightSim.py \
	UCSCScheduler_V3.py \
	fake_apflog.py \
	sim_night.py \
	sim_nights.py \
	x_gaussslit.py \
	getpriority.py \
	gen_dates.py \
	make_sim_files.py \
	makevels.py \
	fillin_phasebins.py \
	Generate_Velocities.py \
	consts.py \
	TESSAPF_assess.py \
	prep_sims.py \
	count_planets.py \
	sim_timerange.py \
	sim_timerange_radvel.py \
	sim_timerange_wbreakpoints.py \
	fit_APF_TESS_radvel.py \
	mcmc_radvel_sims.py \
	make_radvel_catalogs.py \
	timedependent_mcmc_radvelfits.py \
	Radvel_MCMC.py \

FITSRC = fit_TESS_APF.py \
	mcmc_velfits.py \
	timedependent_mcmc_velfits.py

install: install_src install_fit

install_src:
	$(INSTALL) $(SRC) $(destdir)

install_fit:
	$(INSTALL) $(FITSRC) $(fitdir)
