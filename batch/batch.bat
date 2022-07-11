ifort /nologo /O3 /Qip /fpp /extend-source:132 /threads /MT /libs:static /dll /exe:GeoMaterials64.dll ^
../code/String.f90            ^
../code/GlobalConstants.f90   ^
../code/Debug.f90             ^
../code/MathLib.f90           ^
../code/StressStrainLib.f90   ^
../code/YieldSurfaceLib.f90   ^
../code/SolversLib.f90        ^
../code/ElastoPlastic.f90     ^
../code/ViscoPlastic.f90      ^
../code/userAdd.f90           ^
../code/userMod.f90

del *.obj
del *.mod
del *.exp
del *.f90
del *.lib
dir *.dll
