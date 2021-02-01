
import geopandas as gpd
import xarray as xr
import pandas as pd
import atlite
import helper
import logging
logger = logging.getLogger(__name__)

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

cutout = atlite.Cutout(snakemake.config['atlite']['cutout_name'],
                       cutout_dir=snakemake.config['atlite']['cutout_dir'])


clustered_busregions_as_geopd = gpd.read_file(
    snakemake.input.regions_onshore).set_index(
        'name', drop=True)

clustered_busregions = pd.Series(clustered_busregions_as_geopd.geometry, index=clustered_busregions_as_geopd.index)

clean_invalid_geometries(clustered_busregions)

I = cutout.indicatormatrix(clustered_busregions)


items = ["total", "urban", "rural"]

pop = pd.DataFrame(columns=items,
                   index=clustered_busregions.index)


for item in items:
    pop_layout = xr.open_dataarray(snakemake.input['pop_layout_' + item])
    pop[item] = I.dot(pop_layout.stack(spatial=('y', 'x')))

pop.to_csv(snakemake.output.clustered_pop_layout)
