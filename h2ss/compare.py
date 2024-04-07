"""Functions to compare and validate results.

References
----------
.. [#Deane21] Deane, P. (2021) Our Climate Neutral Future: Zero by 50. Wind
    Energy Ireland. Available at:
    https://windenergyireland.com/images/files/our-climate-neutral-future-0by50-final-report.pdf
    (Accessed: 8 February 2024).
.. [#Jannel22] Jannel, H. and Torquet, M. (2022). Conceptual design of salt
    cavern and porous media underground storage site. Hystories deliverable
    D7.1-1. Hystories. Available at:
    https://hystories.eu/wp-content/uploads/2022/05/Hystories_D7.1-1-Conceptual-design-of-salt-cavern-and-porous-media-underground-storage-site.pdf
    (Accessed: 9 October 2023).
"""


def electricity_demand_ie(caverns_df):
    """Compare the total capacity to Ireland's electricity demand in 2050.

    Parameters
    ----------
    cavern_df : geopandas.GeoDataFrame
        Geodataframe of caverns within the zone of interest

    Notes
    -----
    Figures from [#Deane21]_.
    """
    print(
        "Energy capacity as a percentage of Ireland's electricity demand "
        "in 2050:",
        f"{(caverns_df['capacity'].sum() / 1000 / 122 * 100):.2f}–"
        f"{(caverns_df['capacity'].sum() / 1000 / 84 * 100):.2f}%"
    )


def cavern_volumes(
    cavern_df, volume_case=380000, minimum_fraction=0.85
):
    """Verify whether cavern volumes are within recommended ranges.

    Parameters
    ----------
    cavern_df : geopandas.GeoDataFrame
        Geodataframe of caverns within the zone of interest
    volume_case : float
        Cavern volume corresponding to a Hystories Project investment scenario
    minimum_fraction : float
        The fraction of ``volume_case`` that is allowed as the minimum

    Returns
    -------
    geopandas.GeoDataFrame
        Dataframe of available caverns

    Notes
    -----
    See [#Jannel22]_ for the Hystories Project investment scenarios. The
    volume used is the free gas volume of a cavern. The volume should be no
    less than 85% of the case's volume.
    """
    return cavern_df[
        cavern_df["cavern_volume"] >= volume_case * minimum_fraction
    ]
