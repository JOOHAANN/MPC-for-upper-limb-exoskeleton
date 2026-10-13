# Compact raw-sensor benchmark cache

`raw_sensor_trials.mat` contains 25 timestamp-aligned trial structures made by
`experiments/bse_generalization/build_bse_trials.m` from the local raw recordings.
It contains recorded kinematics, causal derivatives, filtered EMG, a fixed EMG
history, reference torque, and trial metadata. It does not contain any fitted
MAP prediction or preprocessing parameter learned from target torque.

`provenance.json` identifies the 75 source `.sto` files by SHA-256 and records
the cache hash. Rebuilding with different raw data requires regenerating that
provenance record; do not treat the old hash record as applying to a new cache.

The cache is an offline research dataset derived from `240528Pilot`. The
reference torque's upstream provenance is not established by this cleanup.
