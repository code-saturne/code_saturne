SET CFDSTUDY_ROOT_DIR=C:\Program Files\Code_Saturne\3.2

SET PYTHONPATH=C:\Program Files\Code_Saturne\3.2\lib\python2.7\site-packages\code_saturne;%PYTHONPATH%
SET PYTHONPATH=C:\Program Files\Code_Saturne\3.2\lib\python2.7\site-packages\salome;%PYTHONPATH%

SET PATH=C:\Program Files\Code_Saturne\3.2\bin;%PATH%
SET PATH=C:\Program Files\Code_Saturne\3.2\lib\salome;%PATH%

CALL C:\SALOME\7.2.0\run_salome.bat --modules=GEOM,SMESH,MED,YACS,PARAVIS,CFDSTUDY
