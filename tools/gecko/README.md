# IMPORTANT INFORMATION FOR SYSTEM ADMINISTRATORS

The tool GECKO is disk-intensive. This means that the algorithm will write a lot to the hard disk drive, which, even though buffered, can affect the performance of other processes that are sharing the same filesystem by starving their access to the disk. 

## Recommended use in concurrent systems

If GECKO is being used concurrently, it is recommended that the input file size is limited in order to avoid overall deterioration of the performance of the system. See details below for an estimate.

## Details

The most disk-consuming step is the writing of "seeds" or "hits". These depend on the length and the similarity of the input sequences. 

- A regular bacterial comparison (1 to 10 MB input sequences) can write around 1 to 5 GB.
- A medium-sized comparison (around 20 or 30 MB input sequences) can write around 10 GB.
- A chromosome comparison (50 to 200 MB input sequences) between human and gorilla can write around 600 GB. Other "less" similar comparison can write around 100 GB or less.

