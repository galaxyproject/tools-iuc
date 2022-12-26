# InterProScan

## Licensed software/data

The standard install of InterProScan does not include the following components, because their license does not allow redistribution:

- SignalP
- TMHMM
- Phobius
- SMART data (THRESHOLD file)

As a Galaxy instance admin, you can install those manually if needed, as long as you (and the Galaxy users) respect the respective licenses.

## Installing licensed components manually

All the following steps assume that:

- You agree with each software license
- InterProScan was installed using Conda in the `IPRSCAN_DIR` directory
- You have run the InterProScan data manager which placed the data files in the `IPRSCAN_DATA_DIR` directory (e.g. `/data/db/data_managers/interproscan/5.59-91.0/`)

Everytime you upgrade InterProScan, you'll need to do the same things: modify the `${IPRSCAN_DIR}/share/InterProScan/interproscan.properties` file + add the SMART file.

### SignalP

Download and install SignalP 4.1 from https://services.healthtech.dtu.dk/service.php?SignalP-4.1

Modify the `${IPRSCAN_DIR}/share/InterProScan/interproscan.properties` file to use your SignalP install (ie the downloaded archive unzipped in `/path/to/signalp/4.1/`):

```
binary.signalp.path=/path/to/signalp/4.1/signalp
signalp.perl.library.dir=/path/to/signalp/4.1/lib/
```

Make sure the SignalP script begins like this:

```
#!/usr/bin/env perl
```

And make sure, you have properly set `$ENV{SIGNALP}` in the SignalP script (around line 13).

### TMHMM

Download and install TMHMM 2.0c from https://services.healthtech.dtu.dk/service.php?TMHMM-2.0

Modify the `${IPRSCAN_DIR}/share/InterProScan/interproscan.properties` file to use your TMHMM install (ie the downloaded archive unzipped in `/path/to/tmhmm/2.0c/`):

```
binary.tmhmm.path=/path/to/tmhmm/2.0c/bin/decodeanhmm.Linux_x86_64
tmhmm.model.path=/path/to/tmhmm/2.0c/lib/TMHMM2.0.model
```

### Phobius

Download Phobius from https://phobius.sbc.su.se/data.html

Modify the `${IPRSCAN_DIR}/share/InterProScan/interproscan.properties` file to use your Phobius install (ie the downloaded archive unzipped in `/pth/to/phobiu/`):

```
binary.phobius.pl.path=/path/to/phobius/phobius.pl
```

Make sure the phobius.pl script begins like this:

```
#!/usr/bin/env perl
```

### SMART data

Download SMART data from https://software.embl-em.de/software/18 (choose SMART 7.1)

Copy the THRESHOLDS file from the archive into `${IPRSCAN_DATA_DIR}/data/smart/7.1/THRESHOLDS`

Modify the `${IPRSCAN_DIR}/share/InterProScan/interproscan.properties` file to use this SMART file:

```
smart.threshold.path=${IPRSCAN_DATA_DIR}/data/smart/7.1/THRESHOLDS
```
