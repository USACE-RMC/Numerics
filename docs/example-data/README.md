# Example Data

This folder contains datasets used in the ***Numerics*** library documentation tutorials. Each dataset is sourced from published references and validated against established statistical software (R, rstan, HEC-SSP).

## Datasets

| File | Description | Records | Source |
|------|-------------|---------|--------|
| `tippecanoe-river-streamflow.csv` | Annual peak streamflow, Tippecanoe River near Delphi, IN (Station 43) | 48 | Rao & Hamed (2000), Table 5.1.1 |
| `white-river-nora-floods.csv` | Annual peak streamflow, White River near Nora, IN | 62 | Rao & Hamed (2000), Table 7.1.2 |
| `white-river-mt-carmel-exceedances.csv` | Threshold exceedances (≥50,000 cfs), White River at Mt. Carmel, IN | 281 | Rao & Hamed (2000), Table 8.3.1 |
| `usgs-01562000-streamflow.csv` | Annual peak streamflow with empirical probabilities, USGS 01562000 | 99 | USGS Bulletin 17C test sites |
| `iris-dataset.csv` | Fisher's Iris flower measurements, 3 species | 150 | Fisher (1936) |
| `statistics-samples.csv` | Two general-purpose samples for statistical analysis | 69 | Validated against R (2024) |
| `streamflow-regression.csv` | Regional watershed characteristics and annual peak flows | 25 | Synthetic data modeled after USGS regional regression relationships |
| `flood-damage-glm.csv` | Flood hydraulics and binary damage outcomes | 35 | Synthetic data for logistic regression modeling |

## References

- Rao, A. R., & Hamed, K. H. (2000). *Flood Frequency Analysis*. CRC Press.
- Fisher, R. A. (1936). The use of multiple measurements in taxonomic problems. *Annals of Eugenics*, 7(2), 179-188.
- U.S. Geological Survey. (2018). Guidelines for Determining Flood Flow Frequency — Bulletin 17C.
- R Core Team (2024). R: A Language and Environment for Statistical Computing.

## Loading Data in C#

```cs
// Example: Load CSV data for analysis
using System.IO;
using System.Linq;

string[] lines = File.ReadAllLines("example-data/tippecanoe-river-streamflow.csv");
double[] data = lines
    .Where(line => !line.StartsWith("#") && !string.IsNullOrWhiteSpace(line))
    .Skip(1) // Skip header
    .Select(line => double.Parse(line.Trim()))
    .ToArray();
```
