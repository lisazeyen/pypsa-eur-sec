
configfile: "config.yaml"

wildcard_constraints:
    lv="[a-z0-9\.]+",
    simpl="[a-zA-Z0-9]*",
    clusters="[0-9]+m?",
    year="[0-9]+m?",
    sectors="[+a-zA-Z0-9]+",
    opts="[-+a-zA-Z0-9]*",
    sector_opts="[-+a-zA-Z0-9]*"



subworkflow pypsaeur:
    workdir: "../pypsa-eur"
    snakefile: "../pypsa-eur/Snakefile"
    configfile: "../pypsa-eur/config.yaml"

rule all:
     input: config['summary_dir'] + '/' + config['run'] + '/{year}/graphs/costs.pdf'


rule solve_all_elec_networks:
    input:
        expand(config['results_dir'] + config['run'] + "/{year}/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc",
               **config['scenario'])

rule test_script:
    input:
        expand("resources/heat_demand_urban_elec_s_{clusters}.nc",
                 **config['scenario'])


rule prepare_sector_networks:
    input:
        expand(config['results_dir'] + config['run'] + "/{year}/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc",
                 **config['scenario'])


rule build_population_layouts:
    input:
        # nuts3_shapes=pypsaeur('resources/nuts3_shapes.geojson'),
        nuts3_shapes='../pypsa-eur/resources/nuts3_shapes.geojson',
        urban_percent="data/urban_percent.csv"
    output:
        pop_layout_total="resources/{year}/pop_layout_total.nc",
        pop_layout_urban="resources/{year}/pop_layout_urban.nc",
        pop_layout_rural="resources/{year}/pop_layout_rural.nc"
    script: "scripts/build_population_layouts.py"


rule build_clustered_population_layouts:
    input:
        pop_layout_total="resources/{year}/pop_layout_total.nc",
        pop_layout_urban="resources/{year}/pop_layout_urban.nc",
        pop_layout_rural="resources/{year}/pop_layout_rural.nc",
 	#regions_onshore='../pypsa-eur/resources/2013/regions_onshore_{network}_s{simpl}_{clusters}.geojson'
        regions_onshore=pypsaeur('resources/{year}/regions_onshore_{network}_s{simpl}_{clusters}.geojson')
    output:
        clustered_pop_layout="resources/{year}/pop_layout_{network}_s{simpl}_{clusters}.csv"
    script: "scripts/build_clustered_population_layouts.py"


rule build_heat_demands:
    input:
        pop_layout_total="resources/{year}/pop_layout_total.nc",
        pop_layout_urban="resources/{year}/pop_layout_urban.nc",
        pop_layout_rural="resources/{year}/pop_layout_rural.nc",
 	#regions_onshore='../pypsa-eur/resources/2013/regions_onshore_{network}_s{simpl}_{clusters}.geojson'
        regions_onshore=pypsaeur("resources/{year}/regions_onshore_{network}_s{simpl}_{clusters}.geojson")
    output:
        heat_demand_urban="resources/{year}/heat_demand_urban_{network}_s{simpl}_{clusters}.nc",
        heat_demand_rural="resources/{year}/heat_demand_rural_{network}_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/{year}/heat_demand_total_{network}_s{simpl}_{clusters}.nc"
    script: "scripts/build_heat_demand.py"

rule build_temperature_profiles:
    input:
        pop_layout_total="resources/{year}/pop_layout_total.nc",
        pop_layout_urban="resources/{year}/pop_layout_urban.nc",
        pop_layout_rural="resources/{year}/pop_layout_rural.nc",
 	#regions_onshore='../pypsa-eur/resources/2013/regions_onshore_{network}_s{simpl}_{clusters}.geojson'
        regions_onshore=pypsaeur("resources/{year}/regions_onshore_{network}_s{simpl}_{clusters}.geojson")
    output:
        temp_soil_total="resources/{year}/temp_soil_total_{network}_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/{year}/temp_soil_rural_{network}_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/{year}/temp_soil_urban_{network}_s{simpl}_{clusters}.nc",
        temp_air_total="resources/{year}/temp_air_total_{network}_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/{year}/temp_air_rural_{network}_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/{year}/temp_air_urban_{network}_s{simpl}_{clusters}.nc"
    script: "scripts/build_temperature_profiles.py"


rule build_cop_profiles:
    input:
        temp_soil_total="resources/{year}/temp_soil_total_{network}_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/{year}/temp_soil_rural_{network}_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/{year}/temp_soil_urban_{network}_s{simpl}_{clusters}.nc",
        temp_air_total="resources/{year}/temp_air_total_{network}_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/{year}/temp_air_rural_{network}_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/{year}/temp_air_urban_{network}_s{simpl}_{clusters}.nc"
    output:
        cop_soil_total="resources/{year}/cop_soil_total_{network}_s{simpl}_{clusters}.nc",
        cop_soil_rural="resources/{year}/cop_soil_rural_{network}_s{simpl}_{clusters}.nc",
        cop_soil_urban="resources/{year}/cop_soil_urban_{network}_s{simpl}_{clusters}.nc",
        cop_air_total="resources/{year}/cop_air_total_{network}_s{simpl}_{clusters}.nc",
        cop_air_rural="resources/{year}/cop_air_rural_{network}_s{simpl}_{clusters}.nc",
        cop_air_urban="resources/{year}/cop_air_urban_{network}_s{simpl}_{clusters}.nc"
    script: "scripts/build_cop_profiles.py"


rule build_solar_thermal_profiles:
    input:
        pop_layout_total="resources/{year}/pop_layout_total.nc",
        pop_layout_urban="resources/{year}/pop_layout_urban.nc",
        pop_layout_rural="resources/{year}/pop_layout_rural.nc",
 	#regions_onshore='../pypsa-eur/resources/2013/regions_onshore_{network}_s{simpl}_{clusters}.geojson'
        regions_onshore=pypsaeur("resources/{year}/regions_onshore_{network}_s{simpl}_{clusters}.geojson")
    output:
        solar_thermal_total="resources/{year}/solar_thermal_total_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_urban="resources/{year}/solar_thermal_urban_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_rural="resources/{year}/solar_thermal_rural_{network}_s{simpl}_{clusters}.nc"
    script: "scripts/build_solar_thermal_profiles.py"



rule build_energy_totals:
    input:
        nuts3_shapes=pypsaeur('resources/nuts3_shapes.geojson'),
	    district_heat_share='data/district_heat_share.csv',
        idee_dir='data/jrc-idees-2015',
        eea_co2="data/eea/UNFCCC_v21.csv",
        swiss="data/switzerland-sfoe/switzerland-new_format.csv"
    output:
        energy_name='data/{year}/energy_totals.csv',
    	co2_name='data/{year}/co2_totals.csv',
    	transport_name='data/{year}/transport_data.csv'
    threads: 1
    resources: mem_mb=10000
    script: 'scripts/build_energy_totals.py'

rule build_biomass_potentials:
    input:
        jrc_potentials="data/biomass/JRC Biomass Potentials.xlsx"
    output:
        biomass_potentials='data/biomass_potentials.csv'
    threads: 1
    resources: mem_mb=1000
    script: 'scripts/build_biomass_potentials.py'


rule build_industry_sector_ratios:
    output:
        industry_sector_ratios="resources/industry_sector_ratios.csv"
    threads: 1
    resources: mem_mb=1000
    script: 'scripts/build_industry_sector_ratios.py'


rule build_industrial_demand_per_country:
    input:
        industry_sector_ratios="resources/industry_sector_ratios.csv"
    output:
        industrial_demand_per_country="resources/industrial_demand_per_country.csv"
    threads: 1
    resources: mem_mb=1000
    script: 'scripts/build_industrial_demand_per_country.py'


rule build_industrial_demand:
    input:
        clustered_pop_layout="resources/pop_layout_{network}_s{simpl}_{clusters}.csv",
        industrial_demand_per_country="resources/industrial_demand_per_country.csv"
    output:
        industrial_demand="resources/industrial_demand_{network}_s{simpl}_{clusters}.csv"
    threads: 1
    resources: mem_mb=1000
    script: 'scripts/build_industrial_demand.py'


if config['sector'].get('retrofitting', True):
	rule build_retro_cost_curves:
	    input:
    		building_stock="data/retro/data_building_stock.csv",
    		u_values_PL="data/retro/u_values_poland.csv",
    		tax_w="data/retro/electricity_taxes_eu.csv",
    		construction_index="data/retro/comparative_level_investment.csv",
    		average_surface="data/retro/average_surface_components.csv",
    		floor_area_missing="data/retro/floor_area_missing.csv",
    		clustered_pop_layout="resources/{year}/pop_layout_{network}_s{simpl}_{clusters}.csv",
    		cost_germany="data/retro/retro_cost_germany.csv",
    		window_assumptions="data/retro/window_assumptions.csv"
	    output:
    		retro_cost="resources/{year}/retro_cost_{network}_s{simpl}_{clusters}.csv",
    		floor_area="resources/{year}/floor_area_{network}_s{simpl}_{clusters}.csv"
	    script: "scripts/build_retro_cost_curves.py"


rule prepare_sector_network:
    input:
        energy_totals_name='data/{year}/energy_totals.csv',
        co2_totals_name='data/{year}/co2_totals.csv',
        transport_name='data/{year}/transport_data.csv',
        biomass_potentials='data/biomass/biomass_potentials.csv',
        biomass_transport='data/biomass/biomass_transport_costs.csv',
        timezone_mappings='data/timezone_mappings.csv',
        heat_profile="data/heat_load_profile_BDEW.csv",
        costs="data/costs/",
        #network='../pypsa-eur/networks/2013/{network}_s{simpl}_{clusters}_lv{lv}_{opts}.nc',
        network=pypsaeur('networks/{year}/{network}_s{simpl}_{clusters}_lv{lv}_{opts}.nc'),
	costs_old="data/costs_old.csv",
	retro_cost_energy = "resources/{year}/retro_cost_{network}_s{simpl}_{clusters}.csv",
        floor_area = "resources/{year}/floor_area_{network}_s{simpl}_{clusters}.csv",
        clustered_pop_layout="resources/{year}/pop_layout_{network}_s{simpl}_{clusters}.csv",
	traffic_data = "data/emobility/",
	h2_cavern = "data/hydrogen_salt_cavern_potentials.csv",
        industrial_demand="resources/industrial_demand_{network}_s{simpl}_{clusters}.csv",
        heat_demand_urban="resources/{year}/heat_demand_urban_{network}_s{simpl}_{clusters}.nc",
        heat_demand_rural="resources/{year}/heat_demand_rural_{network}_s{simpl}_{clusters}.nc",
        heat_demand_total="resources/{year}/heat_demand_total_{network}_s{simpl}_{clusters}.nc",
        temp_soil_total="resources/{year}/temp_soil_total_{network}_s{simpl}_{clusters}.nc",
        temp_soil_rural="resources/{year}/temp_soil_rural_{network}_s{simpl}_{clusters}.nc",
        temp_soil_urban="resources/{year}/temp_soil_urban_{network}_s{simpl}_{clusters}.nc",
        temp_air_total="resources/{year}/temp_air_total_{network}_s{simpl}_{clusters}.nc",
        temp_air_rural="resources/{year}/temp_air_rural_{network}_s{simpl}_{clusters}.nc",
        temp_air_urban="resources/{year}/temp_air_urban_{network}_s{simpl}_{clusters}.nc",
        cop_soil_total="resources/{year}/cop_soil_total_{network}_s{simpl}_{clusters}.nc",
        cop_soil_rural="resources/{year}/cop_soil_rural_{network}_s{simpl}_{clusters}.nc",
        cop_soil_urban="resources/{year}/cop_soil_urban_{network}_s{simpl}_{clusters}.nc",
        cop_air_total="resources/{year}/cop_air_total_{network}_s{simpl}_{clusters}.nc",
        cop_air_rural="resources/{year}/cop_air_rural_{network}_s{simpl}_{clusters}.nc",
        cop_air_urban="resources/{year}/cop_air_urban_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_total="resources/{year}/solar_thermal_total_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_urban="resources/{year}/solar_thermal_urban_{network}_s{simpl}_{clusters}.nc",
        solar_thermal_rural="resources/{year}/solar_thermal_rural_{network}_s{simpl}_{clusters}.nc"
    output:
        network = config['results_dir']  +  config['run'] + '/{year}/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc' #,
	# costs = config['results_dir']  +  config['run'] + '/costs/assumed_costs_{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.csv'
    threads: 1
    resources: mem=1000
    benchmark: "benchmarks/prepare_network/{year}/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}"
    script: "scripts/prepare_sector_network.py"


rule solve_network:
    input:
        network=config['results_dir'] + config['run'] + "/{year}/prenetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc",
        config=config['summary_dir'] + '/' + config['run'] + '/configs/config.yaml'
    output: config['results_dir'] + config['run'] + "/{year}/postnetworks/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc"
    shadow: "shallow"
    log:
        solver="logs/" + config['run'] + "/{year}/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_solver.log",
        python="logs/" + config['run'] + "/{year}/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_python.log",
        memory="logs/" + config['run'] + "/{year}/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}_memory.log"
    benchmark: "benchmarks/solve_network/{year}/{network}_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}"
    threads: 4
    resources: mem=100000 #memory in MB; 40 GB enough for 45+B+I; 100 GB based on RESI usage for 128
    # group: "solve" # with group, threads is ignored https://bitbucket.org/snakemake/snakemake/issues/971/group-job-description-does-not-contain
    script: "scripts/solve_network.py"

rule plot_network:
    input:
        network=config['results_dir'] + config['run'] + "/{year}/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc"
    output:
        map=config['results_dir'] + config['run'] + "/{year}/maps/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}-costs-all.pdf",
	today=config['results_dir'] + config['run'] + "/{year}/maps/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}-today.pdf"
    threads: 2
    resources: mem_mb=10000
    script: "scripts/plot_network.py"



rule copy_config:
    input:
        networks=expand(config['results_dir'] + config['run'] + "/{year}/prenetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc",
               **config['scenario'])
    output:
        config=config['summary_dir'] + '/' + config['run'] + '/configs/config.yaml'
    threads: 1
    resources: mem_mb=1000
    script:
        'scripts/copy_config.py'


rule make_summary:
    input:
        networks=expand(config['results_dir'] + config['run'] + "/{year}/postnetworks/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}.nc",
               **config['scenario']),
        plots=expand(config['results_dir'] + config['run'] + "/{year}/maps/elec_s{simpl}_{clusters}_lv{lv}_{opts}_{sector_opts}-costs-all.pdf",
               **config['scenario'])
        #heat_demand_name='data/heating/daily_heat_demand.h5'
    output:
        nodal_costs=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/nodal_costs.csv',
        nodal_capacities=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/nodal_capacities.csv',
        nodal_cfs=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/nodal_cfs.csv',
        cfs=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/cfs.csv',
        costs=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/costs.csv',
        capacities=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/capacities.csv',
        curtailment=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/curtailment.csv',
        energy=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/energy.csv',
        supply=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/supply.csv',
        supply_energy=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/supply_energy.csv',
        prices=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/prices.csv',
        weighted_prices=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/weighted_prices.csv',
        market_values=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/market_values.csv',
        price_statistics=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/price_statistics.csv',
        metrics=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/metrics.csv'
    threads: 2
    resources: mem_mb=10000
    script:
        'scripts/make_summary.py'


rule plot_summary:
    input:
        costs=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/costs.csv',
        energy=config['summary_dir'] + '/' + config['run'] + '/{year}/csvs/energy.csv'
    output:
        costs=config['summary_dir'] + '/' + config['run'] + '/{year}/graphs/costs.pdf',
        energy=config['summary_dir'] + '/' + config['run'] + '/{year}/graphs/energy.pdf'
    threads: 2
    resources: mem_mb=10000
    script:
        'scripts/plot_summary.py'
