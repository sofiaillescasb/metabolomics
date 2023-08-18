@echo off

set "sourceDir=C:\Users\killescas\OneDrive - San Juan de Dios\Desktop\GitHub\metabolomics"
set "destinationDir=\\hsjdbcn.es\dfsroot\Recursos\metabolismosinaptico\SOFIA\Omics\01.neurodevelopmental_disorders\metabolomics"

robocopy "%sourceDir%" "%destinationDir%" /E /Z /COPY:DAT /R:3 /W:5 /NP /XO /XX /IS /IT


set "sourceDir2=C:\Users\killescas\OneDrive - San Juan de Dios\Desktop\GitHub\lipidomics"
set "destinationDir2=\\hsjdbcn.es\dfsroot\Recursos\metabolismosinaptico\SOFIA\Omics\01.neurodevelopmental_disorders\lipidomics"

robocopy "%sourceDir2%" "%destinationDir2%" /E /Z /COPY:DAT /R:3 /W:5 /NP /XO /XX /IS /IT