
import geopandas as gpd
import atlite
import pandas as pd
import xarray as xr
import scipy as sp
import helper

if 'snakemake' not in globals():
    from vresutils import Dict
    import yaml
    snakemake = Dict()
    with open('config.yaml') as f:
        snakemake.config = yaml.load(f)
    snakemake.input = Dict()
    snakemake.output = Dict()

# adjust snapshots to energy year
snakemake.config["snapshots"] = {'start': '{}-01-01'.format(snakemake.wildcards["year"]),
                                 'end': '{}-01-01'.format(str(int(snakemake.wildcards["year"])+1)),
                                 'closed': 'left'}
snakemake.config['atlite']['cutout_name'] = 'europe-{}'.format(snakemake.wildcards.year)

time = pd.date_range(freq='m', **snakemake.config['snapshots'])
params = dict(years=slice(*time.year[[0, -1]]),
              months=slice(*time.month[[0, -1]]))

cutout = atlite.Cutout(snakemake.config['atlite']['cutout_name'],
                       cutout_dir=snakemake.config['atlite']['cutout_dir'],
                       **params)

clustered_busregions_as_geopd = gpd.read_file(
    snakemake.input.regions_onshore).set_index(
        'name', drop=True)

clustered_busregions = pd.Series(
    clustered_busregions_as_geopd.geometry,
    index=clustered_busregions_as_geopd.index)

helper.clean_invalid_geometries(clustered_busregions)

I = cutout.indicatormatrix(clustered_busregions)


for item in ["rural", "urban", "total"]:

    pop_layout = xr.open_dataarray(snakemake.input['pop_layout_'+item])

    M = I.T.dot(sp.diag(I.dot(pop_layout.stack(spatial=('y', 'x')))))

    heat_demand = cutout.heat_demand(matrix=M.T,index=clustered_busregions.index)

    heat_demand = cutout.heat_demand(
        matrix=M.T, index=clustered_busregions.index)

    heat_demand.to_netcdf(snakemake.output["heat_demand_" + item])
