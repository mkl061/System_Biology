For running the R file you need these, in following order:

R (core): https://www.r-project.org/

RStudio: https://www.rstudio.com/products/rstudio/download/

RTools: https://cran.r-project.org/bin/windows/Rtools/


OBS! The following isd based on Windows, so use the Mac-version of Ctrl (i.e Cmd)

1. Download the latest versions of the Excel sheets "nodes" and "edges" as CSV, naming them "nodes" and "edges".
2. Open Cytoscape
3. Now you can open the R file in RStudio.
4. Then click Ctrl+Shift+H, and locate the folder the to CSV files (nodes.csv, edges.csv)
5. Lastly, select the entire script by clicking Ctrl+A, followed by Ctrl+Enter
6. Everything should be updated in Cytoscape when the script is done loading!


OBS! The script may take a long time to run the first time because it needs to download packages. 
If nothing happens in Cytoscape the most likely reason is that the CSV file has a separator other than
";" (which is the Norwegian default). The other separator is ",". So then you need to change the script.
It should be clearly marked where.


The packages used for this project are:
- Tidyverse (mainly; dplyr, tidyr, stringr)
- BiocManager (i.e. BioConductor)
- biomaRt (i.e. Ensembl BioMart)