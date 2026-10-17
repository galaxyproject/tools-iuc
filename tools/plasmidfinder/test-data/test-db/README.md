The test database contains three reference sequences. Rebuild its KMA indexes
from this directory with KMA 1.6.18:

```sh
for database in Inc18 Rep3 Rep_trans; do
    kma index -k 16 -m 5 -i "${database}.fsa" -o "${database}"
done
```
