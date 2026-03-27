import cfr
import xarray as xr
import pandas as pd
import numpy as np
import sys

pdb = cfr.ProxyDatabase().fetch('PAGES2kv2')
df = pdb.to_df()


# Example: df has columns
# ['pid', 'lat', 'lon', 'elev', 'ptype', 'time', 'value']

# Step 1: explode list-columns so each proxy-year is one row
long_df = df[['pid', 'lat', 'lon', 'elev', 'ptype', 'time', 'value']].explode(
    ['time', 'value'],
    ignore_index=True
)

# Step 2: make sure types are numeric
long_df['time'] = pd.to_numeric(long_df['time'], errors='coerce')
long_df['value'] = pd.to_numeric(long_df['value'], errors='coerce')
long_df['year'] = np.floor(long_df['time']).astype('Int64')

# Step 3: create full year index, e.g. 1 to 2000
full_years = pd.Index(range(1, 2001), name='year')

# Step 4: pivot to wide matrix (and account for multiple observations per year by taking the mean)
proxy_matrix = (
    long_df
    .assign(year=np.floor(long_df['time']).astype('Int64'))
    .pivot_table(index='year', columns='pid', values='value', aggfunc='mean')
    .reindex(range(1, 2001))
)

# Step 5: metadata table
proxy_meta = (
    df[['pid', 'lat', 'lon', 'elev', 'ptype']]
    .drop_duplicates('pid')
    .set_index('pid')
    .loc[proxy_matrix.columns]
)

# Step 6: combine metadata and data into one DataFrame (optional)
meta_rows = pd.DataFrame(
    [proxy_meta.loc[proxy_matrix.columns, 'lat'],
     proxy_meta.loc[proxy_matrix.columns, 'lon'],
     proxy_meta.loc[proxy_matrix.columns, 'elev'],
     proxy_meta.loc[proxy_matrix.columns, 'ptype']],
    index=['lat', 'lon', 'elev', 'ptype']
)

combined = pd.concat([meta_rows, proxy_matrix])
print(combined)
proxy_matrix.to_csv('data/PAGES2K_proxy_matrix_1-2000.csv', index=True)
proxy_meta.to_csv('data/PAGES2K_proxy_metadata_1-2000.csv', index=True)
combined.to_csv('data/PAGES2K_proxy_combined_data_1-2000.csv', index=True)