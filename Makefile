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
	sim_timerange.py

FITSRC = fit_TESS_APF.py

install: install_src install_fit

install_src:
	$(INSTALL) $(SRC) $(destdir)

install_fit:
	$(INSTALL) $(FITSRC) $(fitdir)
