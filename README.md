# frib_hrs_drift_chambers
Garfield++, Gmsh, and Elmer simulations of drift chambers for the FRIB HRS drift chambers.

## Elmer field maps

The solved Elmer mesh and field results for each cage are stored as zip files:

* `StripCage/strip_cage.zip` is included in this repository.
* `field_cage.zip` is over GitHub's file size limit, so it is attached to the [fieldcage release]).
* 
Before running either simulation, unzip `strip_cage.zip` inside `StripCage/`, and download `field_cage.zip` from the release and unzip it inside `WFCDC/`, so the mesh and result files sit where each program expects them.
