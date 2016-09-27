download=wget

all: rignot bedmap2 searise arthern

.PHONY: le-brocq searise rignot rignot-hires bedmap2 bedmap2-rignot arthern


searise: Antarctica_5km_withshelves_v0.75.nc

Antarctica_5km_withshelves_v0.75.nc:
	$(download) http://websrv.cs.umt.edu/isis/images/c/cc/Antarctica_5km_withshelves_v0.75.nc


le-brocq: Antarctica_5km_dev1.0.nc

Antarctica_5km_dev1.0.nc:
	$(download) http://websrv.cs.umt.edu/isis/images/4/4d/Antarctica_5km_dev1.0.nc


rignot: antarctica_ice_velocity_900m.nc

rignot-hires: antarctica_ice_velocity_450m.nc

antarctica_ice_velocity_900m.nc:
	$(download) ftp://n5eil01u.ecs.nsidc.org/SAN/MEASURES/NSIDC-0484.001/1996.01.01/antarctica_ice_velocity_900m.nc

antarctica_ice_velocity_450m.nc:
	$(download) ftp://n5eil01u.ecs.nsidc.org/SAN/MEASURES/NSIDC-0484.001/1996.01.01/antarctica_ice_velocity_450m.nc


bedmap2: bedmap2_bin.zip

bedmap2_bin.zip:
	$(download) https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_bin.zip
	unzip $@
	$(download) https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_readme.txt
	$(download) https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_Updates.rtf
	$(download) https://secure.antarctica.ac.uk/data/bedmap2/resources/matlab/load_bedmap2.m

bedmap2-rignot: bin.zip

bin.zip:
	$(download) https://secure.antarctica.ac.uk/data/bedmap2/resources/Rignot_velocity/bin.zip
	unzip $@
	$(download) https://secure.antarctica.ac.uk/data/bedmap2/resources/Rignot_velocity/Rignot_velocity_bedmap2_grid_readme.rtf

arthern: Arthern_accumulation_bin.zip

Arthern_accumulation_bin.zip:
	$(download) https://secure.antarctica.ac.uk/data/bedmap2/resources/Arthern_accumulation/Arthern_accumulation_bedmap2_grid_readme.rtf
	$(download) https://secure.antarctica.ac.uk/data/bedmap2/resources/Arthern_accumulation/Arthern_accumulation_bin.zip
	unzip $@
