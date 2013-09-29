SET CSVERS=3.2-alpha

SET CFDSTUDY_ROOT_DIR=C:\Program Files\Code_Saturne\%CSVERS%

SET PYTHONPATH=C:\Program Files\Code_Saturne\%CSVERS%\lib\python2.7\site-packages\code_saturne;%PYTHONPATH%
SET PYTHONPATH=C:\Program Files\Code_Saturne\%CSVERS%\lib\python2.7\site-packages\salome;%PYTHONPATH%

SET PATH=C:\Program Files\Code_Saturne\%CSVERS%\bin;%PATH%
SET PATH=C:\Program Files\Code_Saturne\%CSVERS%\lib\salome;%PATH%

CALL C:\SALOME\7.2.0\run_salome.bat --modules=GEOM,SMESH,MED,YACS,PARAVIS,CFDSTUDY
