The test database contains three reference sequences. Rebuild its KMA indexes
from this directory with KMA 1.6.18:

```sh
for database in Inc18 Rep3 Rep_trans; do
    kma index -i "${database}.fsa" -o "${database}"
done
```

Each index must contain only the sequences from its corresponding FASTA file.
