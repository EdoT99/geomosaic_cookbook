# This is a cookbook for downstream analysis from  Geomosaic pipeline [Corso D. et al, 2024]
# 🧬 Geomosaic Downstream Analysis Notebook
---

### 👤 Authorship
* **Edoardo Taccaliti** ([edoardotaccaliti@gmail.com](mailto:edoardotaccaliti@gmail.com))
* **Davide Corso** ([davidecrs92@gmail.com](mailto:davidecrs92@gmail.com))

### 📜 Copyright & Licensing
**Copyright (c) 2026 Edoardo Taccaliti & Davide Corso**


### Set up
Create the environment
```bash 
mamba env create -f env.yaml 

```
Activate the environemnt:
```bash 
mamba env activate geomosaic_analysis

```
IMPORTANT: if you use RStudio you need to [launch it from the terminal](https://stackoverflow.com/questions/38534383/how-to-set-up-conda-installed-r-for-use-with-rstudio/62737170#62737170) with the activated conda environment:

 ```bash
mamba env activate geomosaic_analysis
rstudio
```

 
This notebook is licensed under the **MIT License**.  
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files...

> **Note:** Please cite the original Geomosaic paper if using this pipeline for published research.

Corso D., et al., "Geomosaic: a flexible metagenomic pipeline combining biological and geochemical data to outline biosphere and geosphere interactions."

