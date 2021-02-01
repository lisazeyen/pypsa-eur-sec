
import geopandas as gpd
import atlite
import pandas as pd
import xarray as xr
import scipy as sp
import helper
import logging
logger = logging.getLogger(__name__)

if 'snakemake' not in globals():
    from vresutils import Dict
    import yaml
    snakemake = Dict()
    with open('config.yaml') as f:
        snakemake.config = yaml.load(f)
    snakemake.input = Dict()
    snakemake.output = Dict()


def clean_invalid_geometries(geometries):
    """Fix self-touching or self-crossing polygons; these seem to appear
    due to numerical problems from writing and reading, since the geometries
    are valid before being written in pypsa-eur/scripts/cluster_network.py"""
    for i,p in geometries.items():
        if not p.is_valid:
            logger.warning(f'Clustered region {i} had an invalid geometry, fixing using zero buffer.')
            geometries[i] = p.buffer(0)

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

clean_invalid_geometries(clustered_busregions)

I = cutout.indicatormatrix(clustered_busregions)


for item in ["total", "rural", "urban"]:

    pop_layout = xr.open_dataarray(snakemake.input['pop_layout_'+item])

    M = I.T.dot(sp.diag(I.dot(pop_layout.stack(spatial=('y', 'x')))))
    nonzero_sum = M.sum(axis=0, keepdims=True)
    nonzero_sum[nonzero_sum == 0.] = 1.
    M_tilde = M/nonzero_sum

    temp_air = cutout.temperature(matrix=M_tilde.T,index=clustered_busregions.index)

    temp_air.to_netcdf(snakemake.output["temp_air_" + item])

    temp_soil = cutout.soil_temperature(
        matrix=M_tilde.T, index=clustered_busregions.index)

    temp_soil.to_netcdf(snakemake.output["temp_soil_" + item])
