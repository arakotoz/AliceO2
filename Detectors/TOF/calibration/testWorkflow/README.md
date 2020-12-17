<!-- doxy
\page refDetectorsTOFtestWorkflow testWorkflow
/doxy -->

# TOF calibration workflows

## DCS DP processing:

Local example workflow with local CCDB (running on port 8080) :

This will read the list of DPs to be associated to TOF from CCDB (remove
`--use-ccdb-to-configure` if you don't want this, but use hardcoded
aliases. You can specify the path of CCDB also with `--ccdb-path`.
YOu can also specify to run in verbose mode (`--use-verbose-mode`)

```shell
o2-dcs-sim-workflow --max-timeframes 3 --delta-fraction 0.5 -b |
o2-calibration-tof-dcs-workflow --use-ccdb-to-configure -b |
o2-calibration-ccdb-populator-workflow --ccdb-path="http://localhost:8080" -b
```
To populate locally a DCS entry for the configuration, run the:

`O2/Detectors/TOF/calibration/macros/makeCCDBEntryForDCS.C`

macro.


## LHC phase:

This will process the LHC phase simulated by the Generator workflow

```shell
LHC phase
o2-calibration-data-generator-workflow --lanes 10 --mean-latency 100000 --max-timeframes 500 |
o2-calibration-lhc-clockphase-workflow --tf-per-slot 20 |
o2-calibration-ccdb-populator-workflow --ccdb-path localhost:8080
```

## TOF channel calibration:

To obtain the TOF channel offsets (at end of processing). Input from simulation, but is should work if attached to reco+calib flow

* simulating reading from ccdb, and using it, with "-b", in "test" mode --> to use this, we need an appropriate CCDB object in the CCDB

```shell
o2-calibration-data-generator-workflow --lanes 10 --mean-latency 100000 --max-timeframes 500 --do-TOF-channel-calib --do-TOF-channel-calib-in-test-mode -b |
o2-testworkflows-tof-dummy-ccdb -b |
o2-calibration-tof-channel-calib-workflow --min-entries 50 --do-TOF-channel-calib-in-test-mode --use-ccdb -b |
o2-calibration-ccdb-populator-workflow --ccdb-path localhost:8080 -b
```

* simulating reading from ccdb, but not using it, with "-b", in "test" mode --> to use this, we need an appropriate CCDB object in the CCDB

```
o2-calibration-data-generator-workflow --lanes 10 --mean-latency 100000 --max-timeframes 500 --do-TOF-channel-calib --do-TOF-channel-calib-in-test-mode -b |
o2-testworkflows-tof-dummy-ccdb -b
| o2-calibration-tof-channel-calib-workflow --min-entries 50 --do-TOF-channel-calib-in-test-mode -b
| o2-calibration-ccdb-populator-workflow --ccdb-path localhost:8080 -b
```

* Using the workflow that has both LHCclockPhase and TOFChannelCalib; for now I can enable only one, or the CCDB populator will not work

```shell
o2-calibration-data-generator-workflow --lanes 10 --mean-latency 100000 --max-timeframes 500 --do-TOF-channel-calib --do-TOF-channel-calib-in-test-mode -b |
o2-calibration-tof-dummy-ccdb-for-calib -b |
o2-calibration-tof-calib-workflow --do-channel-offset --min-entries 50 --do-TOF-channel-calib-in-test-mode -b |
o2-calibration-ccdb-populator-workflow --ccdb-path localhost:8080 -b
```
* same as above, enabling CCDB

```shell
o2-calibration-data-generator-workflow --lanes 10 --mean-latency 100000 --max-timeframes 500 --do-TOF-channel-calib --do-TOF-channel-calib-in-test-mode -b |
o2-calibration-tof-dummy-ccdb-for-calib -b |
o2-calibration-tof-calib-workflow --do-channel-offset --use-ccdb --min-entries 50 --do-TOF-channel-calib-in-test-mode -b |
o2-calibration-ccdb-populator-workflow --ccdb-path localhost:8080 -b
```

## TimeSlewing:

For Time Slewing. Will save the Time Slewing information in files when a certain condition is reached. A post-processing
should then take care of extracting the CCDB

* test mode:

``` shell
o2-calibration-data-generator-workflow --lanes 10 --mean-latency 100000 --max-timeframes 10 -b |
o2-calibration-tof-collect-calib-workflow --tf-sending-policy --running-in-test-mode -b

* non-test but simplified (using option "is-max-number-hits-to-fill-tree-absolute"):

```shell
o2-calibration-data-generator-workflow --lanes 10 --mean-latency 100000 --max-timeframes 10 -b |
o2-calibration-tof-collect-calib-workflow --max-number-hits-to-fill-tree 300 --is-max-number-hits-to-fill-tree-absolute -b
```