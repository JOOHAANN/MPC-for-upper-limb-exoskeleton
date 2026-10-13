# Raw recordings

The local `240528Pilot/` directory contains the supplied acquisition files.
Large recordings are excluded from Git. For a new checkout, place this directory
here, or set `EXO_DATA_ROOT` to a directory containing `240528Pilot`.

```matlab
setenv('EXO_DATA_ROOT', 'D:/datasets/exoskeleton');
```

The path resolver supports both `240528Pilot/CalibrationData/<trial>` and the
supplied `240528Pilot/<trial>` layout. It reports missing files explicitly.
Sample ranges and signal-channel selections remain experiment-specific.
