"""Functions to compare and validate results.

References
----------
.. [#Deane21] Deane, P. (2021) Our Climate Neutral Future: Zero by 50. Wind
    Energy Ireland. Available at:
    https://windenergyireland.com/images/files/our-climate-neutral-future-0by50-final-report.pdf
    (Accessed: 8 February 2024).
"""

def electricity_demand_ie(caverns_df):
    """Compare the total capacity to Ireland's electricity demand in 2050.

    Notes
    -----
    Figures from [#Deane21]_.
    """
    print(
        f"""
            "Energy capacity as a percentage of Ireland's electricity demand
            in 2050:
            {(caverns_df['capacity'].sum() / 1000 / 122 * 100):.2f}–
            {(caverns_df['capacity'].sum() / 1000 / 84 * 100):.2f}%,
        """
    )
